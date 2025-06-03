import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import meshio

def create_mesh(width=20, height=50, mesh_size=0.5, hole_centers=None, hole_radii=None):
    """Create a mesh with circular holes using Gmsh."""
    import gmsh

    gmsh.initialize()
    gmsh.model.add("viscoelastic_plate")

    # Create the rectangle
    rect = gmsh.model.occ.addRectangle(0, 0, 0, width, height)
    
    # Add holes if specified
    holes = []
    if hole_centers is not None and hole_radii is not None:
        for (x, y), r in zip(hole_centers, hole_radii):
            hole = gmsh.model.occ.addDisk(x, y, 0, r, r)
            holes.append(hole)
    
    # Cut holes from the rectangle
    if holes:
        gmsh.model.occ.cut([(2, rect)], [(2, h) for h in holes])
    
    gmsh.model.occ.synchronize()
    
    # Define physical groups for boundary conditions
    bottom_edge = gmsh.model.getBoundary([(2, rect)], recursive=True)
    top_edge = gmsh.model.getBoundary([(2, rect)], recursive=True)
    
    # Set mesh size
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    
    # Generate mesh
    gmsh.model.mesh.generate(2)
    
    # Save mesh
    gmsh.write("mesh.msh")
    
    # Get mesh data
    nodes = np.array(gmsh.model.mesh.getNodes()[1])
    elements = np.array(gmsh.model.mesh.getElements(2)[2][0])
    
    gmsh.finalize()
    
    return nodes, elements

class ViscoelasticFEM:
    def __init__(self, nodes, elements, prony_series, dt=0.1, T=100):
        """
        Initialize viscoelastic FEM solver
        
        Parameters:
        -----------
        nodes : array_like
            Node coordinates [n_nodes, 2]
        elements : array_like
            Element connectivity [n_elements, 3]
        prony_series : dict
            Prony series coefficients with keys 'g_i', 'k_i', and 'tau_i'
        dt : float
            Time step
        T : float
            Total simulation time
        """
        self.nodes = nodes
        self.elements = elements
        self.n_nodes = len(nodes)
        self.n_elements = len(elements)
        self.n_dofs = 2 * self.n_nodes
        
        # Prony series coefficients
        self.g_inf = 1.0 - sum(prony_series['g_i'])
        self.k_inf = 1.0 - sum(prony_series['k_i'])
        self.g_i = prony_series['g_i']
        self.k_i = prony_series['k_i']
        self.tau_i = prony_series['tau_i']
        self.n_terms = len(self.g_i)
        
        # Material constants
        self.E = 1000.0  # Young's modulus
        self.nu = 0.3    # Poisson's ratio
        
        # Convert to Lamé parameters for plane stress
        self.G0 = self.E / (2 * (1 + self.nu))
        self.K0 = self.E / (3 * (1 - 2 * self.nu))
        
        # Time discretization
        self.dt = dt
        self.T = T
        self.n_steps = int(T / dt) + 1
        self.time = np.linspace(0, T, self.n_steps)
        
        # Initialize history variables
        self.e_v = np.zeros((self.n_elements, 3, self.n_terms))  # Viscous deviatoric strain
        self.theta_v = np.zeros((self.n_elements, 1, self.n_terms))  # Viscous volumetric strain
        
        # Initialize displacement storage
        self.u = np.zeros((self.n_steps, self.n_dofs))
    
    def assemble_system(self, step):
        """Assemble the global stiffness matrix and force vector."""
        K = np.zeros((self.n_dofs, self.n_dofs))
        F = np.zeros(self.n_dofs)
        
        # Define m vector for plane stress
        m = np.array([1, 1, 0])
        
        # Define D_mu matrix for plane stress
        D_mu = np.diag([2, 2, 1])
        
        # Define projectors
        I_dev = np.eye(3) - 1/2 * np.outer(m, m)
        
        for el in range(self.n_elements):
            # Get element nodes
            nodes_idx = self.elements[el]
            node_coords = self.nodes[nodes_idx]
            
            # Element DOFs
            dofs = np.array([[2*idx, 2*idx+1] for idx in nodes_idx]).flatten()
            
            # Calculate element area and shape function derivatives
            x = node_coords[:, 0]
            y = node_coords[:, 1]
            
            # Area of the triangle
            area = 0.5 * abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))
            
            # Shape function derivatives
            b = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]]) / (2 * area)
            c = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]]) / (2 * area)
            
            # B matrix (strain-displacement)
            B = np.zeros((3, 6))
            B[0, 0::2] = b
            B[1, 1::2] = c
            B[2, 0::2] = c
            B[2, 1::2] = b
            
            # Deviatoric B matrix
            B_D = I_dev @ B
            
            # Volumetric operator
            b_vol = m.T @ B
            
            # Time-dependent stiffness components
            G_t = self.G0 * (self.g_inf + sum(self.g_i[i] * (1 - self.dt/(2*self.tau_i[i] + self.dt)) 
                                           for i in range(self.n_terms)))
            K_t = self.K0 * (self.k_inf + sum(self.k_i[i] * (1 - self.dt/(2*self.tau_i[i] + self.dt)) 
                                           for i in range(self.n_terms)))
            
            # Element stiffness matrix
            K_el = (G_t * B.T @ D_mu @ B_D + K_t * b_vol.T @ b_vol) * area
            
            # Element history force vector
            F_hist = np.zeros(6)
            if step > 0:
                for i in range(self.n_terms):
                    # Update history variables
                    beta_D = (self.dt/(2*self.tau_i[i] + self.dt) * 
                             (B_D @ self.u[step-1, dofs] - self.e_v[el, :, i]) + 
                             (2*self.tau_i[i])/(2*self.tau_i[i] + self.dt) * self.e_v[el, :, i])
                    
                    beta_V = (self.dt/(2*self.tau_i[i] + self.dt) * 
                             (b_vol @ self.u[step-1, dofs] - self.theta_v[el, 0, i]) + 
                             (2*self.tau_i[i])/(2*self.tau_i[i] + self.dt) * self.theta_v[el, 0, i])
                    
                    # Save updated history variables
                    self.e_v[el, :, i] = (self.dt/(2*self.tau_i[i] + self.dt)) * B_D @ self.u[step, dofs] + beta_D
                    self.theta_v[el, 0, i] = (self.dt/(2*self.tau_i[i] + self.dt)) * b_vol @ self.u[step, dofs] + beta_V
                    
                    # Add to history force vector
                    F_hist += (self.G0 * self.g_i[i] * B.T @ D_mu @ beta_D + 
                              self.K0 * self.k_i[i] * b_vol.T * beta_V)
            
            # Assemble into global system
            for i in range(6):
                F[dofs[i]] += F_hist[i]
                for j in range(6):
                    K[dofs[i], dofs[j]] += K_el[i, j]
        
        return K, F
    
    def apply_boundary_conditions(self, K, F, step):
        """Apply boundary conditions."""
        # Find nodes on bottom and top edges
        bottom_nodes = np.where(np.isclose(self.nodes[:, 1], 0))[0]
        top_nodes = np.where(np.isclose(self.nodes[:, 1], 50))[0]
        
        # Fix bottom edge in Y direction
        for node in bottom_nodes:
            dof = 2 * node + 1  # Y-component
            K[dof, :] = 0
            K[:, dof] = 0
            K[dof, dof] = 1
            F[dof] = 0
        
        # Apply force on top edge (200 N)
        force_per_node = 200 / len(top_nodes)
        for node in top_nodes:
            dof = 2 * node + 1  # Y-component
            F[dof] += force_per_node
        
        return K, F
    
    def solve_time_step(self, step):
        """Solve the system for a single time step."""
        K, F = self.assemble_system(step)
        K, F = self.apply_boundary_conditions(K, F, step)
        
        # Convert to sparse matrix for efficiency
        K_sparse = csr_matrix(K)
        
        # Solve the system
        self.u[step, :] = spsolve(K_sparse, F)
        
        return self.u[step, :]
    
    def solve(self):
        """Solve the viscoelastic problem for all time steps."""
        for step in range(self.n_steps):
            print(f"Solving step {step+1}/{self.n_steps}")
            self.solve_time_step(step)
        
        return self.u
    
    def plot_deformation(self, step):
        """Plot the deformed mesh at a given time step."""
        u_step = self.u[step, :]
        
        # Extract x and y displacements
        u_x = u_step[0::2]
        u_y = u_step[1::2]
        
        # Scale factor for visualization
        scale = 10.0
        
        # Create figure
        plt.figure(figsize=(10, 8))
        
        # Plot undeformed mesh
        for el in range(self.n_elements):
            nodes_idx = self.elements[el]
            x = self.nodes[nodes_idx, 0]
            y = self.nodes[nodes_idx, 1]
            plt.plot(np.append(x, x[0]), np.append(y, y[0]), 'k-', alpha=0.3)
        
        # Plot deformed mesh
        for el in range(self.n_elements):
            nodes_idx = self.elements[el]
            x = self.nodes[nodes_idx, 0] + scale * u_x[nodes_idx]
            y = self.nodes[nodes_idx, 1] + scale * u_y[nodes_idx]
            plt.plot(np.append(x, x[0]), np.append(y, y[0]), 'r-')
        
        plt.axis('equal')
        plt.title(f'Deformation at t = {self.time[step]:.2f} s')
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.grid(True)
        plt.show()
    
    def plot_displacement_history(self):
        """Plot the displacement history at the top center of the specimen."""
        # Find a node near the top center
        center_x = 10
        top_y = 50
        distances = np.sqrt((self.nodes[:, 0] - center_x)**2 + (self.nodes[:, 1] - top_y)**2)
        top_center_node = np.argmin(distances)
        
        # Get displacement history
        top_center_y_dof = 2 * top_center_node + 1
        y_displacement = self.u[:, top_center_y_dof]
        
        # Plot
        plt.figure(figsize=(10, 6))
        plt.plot(self.time, y_displacement)
        plt.title('Displacement History at Top Center')
        plt.xlabel('Time (s)')
        plt.ylabel('Displacement (mm)')
        plt.grid(True)
        plt.show()

def main():
    # Define Prony series coefficients from the table
    prony_series = {
        'g_i': [0.1, 0.08806, 0.07907],
        'k_i': [0.1, 0.08806, 0.07907],
        'tau_i': [1.847, 19.37, 139.9]
    }
    
    # Create mesh with three holes
    width = 20
    height = 50
    hole_centers = [(10, 15), (5, 25), (15, 35)]
    hole_radii = [3, 2, 2.5]
    
    # Try to use gmsh, fall back to simple mesh if not available
    try:
        nodes, elements = create_mesh(width, height, mesh_size=0.5, 
                                      hole_centers=hole_centers, 
                                      hole_radii=hole_radii)
    except ImportError:
        print("Gmsh not available, using simple mesh without holes")
        # Create a simple mesh without holes
        nx, ny = 20, 50
        x = np.linspace(0, width, nx)
        y = np.linspace(0, height, ny)
        X, Y = np.meshgrid(x, y)
        nodes = np.column_stack([X.flatten(), Y.flatten()])
        
        # Create triangular elements
        elements = []
        for j in range(ny-1):
            for i in range(nx-1):
                idx = j * nx + i
                elements.append([idx, idx+1, idx+nx])
                elements.append([idx+1, idx+nx+1, idx+nx])
        elements = np.array(elements)
    
    # Create and solve the FEM problem
    fem = ViscoelasticFEM(nodes, elements, prony_series, dt=1.0, T=500)
    fem.solve()
    
    # Plot results
    fem.plot_deformation(step=50)
    fem.plot_deformation(step=200)
    fem.plot_deformation(step=400)
    fem.plot_displacement_history()

if __name__ == "__main__":
    main()