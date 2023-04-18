// Definition of the Python embedding C API for aims. The implementation is in
// python_interface.py.

struct GridPoint_t {
    double *coords;
    int index_atom, index_radial, index_angular;
};

struct BatchOfPoints_t {
    int size;
    struct GridPoint_t *points;
};

struct AimsContext_t {
    _Bool real_eigenvectors;
    char *elements;
    int n_atoms, n_full_points, n_spin, n_my_batches, n_k_points,
        n_k_points_task, n_basis, n_states, n_hamiltonian_matrix_size;
    int *species;
    double *coords, *rho, *rho_gradient, *kinetic_density,
           *partition_tab, *hirshfeld_volume,
           *hirshfeldw, *freeintegral,
           *KS_eigenvector, *KS_eigenvector_complex, *KS_eigenvalue,
           *occ_numbers, *hamiltonian;
    struct BatchOfPoints_t *batches;
};

extern int call_python_cffi(struct AimsContext_t *, char *, char *, int);
