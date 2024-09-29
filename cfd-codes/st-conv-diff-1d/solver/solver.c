/* Solver for stationary convection-diffusion in 1 dimension 
 * Version: 0.1
 * gcc -o solver .\solver.c
 * .\solver.exe | Out-File log
*/

// Includes
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Requires
#define MESH_FILE "..\\mesher\\meshPoints.txt"
#define SETTINGS_FILE "solverControls.txt"
#define RESULTS_FILE "results.txt"
#define MAXLINE 20
#define U_INIT 1.0
#define PHI_INIT 0.0
#define PHI_LEFT 0.0
#define PHI_RIGHT 1.0
#define RHO 1.0
#define GAMMA 0.02
#define CDS_name "CDS"
#define UDS_name "UDS"
#define CDS 0
#define UDS 1

// Prototypes
void read_mesh(void);
void read_settings(void);
void initialise(void);
void construct_matrices(void);
void solve_tdma(void);
void write_results(void);
void free_mem(void);
void print_matrix(char *s, float *M, int dim_size, int n_dims);
void fprint_matrix(FILE *out, char *s, float *M, int dim_size, int n_dims);
float A_w(int i);
float A_e(int i);
float min(float a, float b);

// Structures
struct solverSettings
{
    float u_init;
    float phi_init;
    float phi_left;
    float phi_right;
    float rho;
    float gamma;
    char scheme_name[MAXLINE];
    int scheme_int;
};

// Global variables
int n_points;
int n_pts_mtrx;
struct solverSettings settings;
float *mesh;
float *u;
float *phi;
float *A;
float *Q;

// Does stuff
int main(void)
{
    printf("Solving a stationary convection-diffusion problem in 1 dimension ...\n");

    // Read mesh
    read_mesh();

    // Read settings
    read_settings();

    // Initialise the solution
    initialise();

    // Construct matrices
    construct_matrices();

    // Manually testing the TDMA solver
    /*
    A[0] = 20;
    A[1] = -5;
    Q[0] = 1100;

    for (int i = 1; i < 4; i++)
    {
        A[i * (n_pts_mtrx + 1) - 1] = -5;
        A[i * (n_pts_mtrx + 1)] = 15;
        A[i * (n_pts_mtrx + 1) + 1] = -5;

        Q[i] = 100;
    }

    A[n_pts_mtrx * n_pts_mtrx - 2] = -5;
    A[n_pts_mtrx * n_pts_mtrx - 1] = 10;
    Q[n_pts_mtrx - 1] = 100;
    */
    /*
    A[0] = 4;
    A[1] = 8;
    Q[0] = 8;

    A[4] = 8;
    A[5] = 18;
    A[6] = 2;
    Q[1] = 18;

    A[9] = 2;
    A[10] = 5;
    A[11] = 1.5;
    Q[2] = 0.5;

    A[14] = 1.5;
    A[15] = 1.75;
    Q[3] = -1.75;

    print_matrix("A:\n", A, n_pts_mtrx, 2);
    print_matrix("Q:\n", Q, n_pts_mtrx, 1);
    */

    // Solve
    solve_tdma();

    // Report
    write_results();
    
    // Cleanup
    free_mem();

    printf("Done!\n");

    return 0;
}

// Reads mesh from MESH_FILE
void read_mesh(void)
{
    printf("Reading mesh file ...\n");

    // Open mesh file
    FILE *fp_mesh = fopen(MESH_FILE, "r");
    if (fp_mesh == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", MESH_FILE);
        exit(2);
    }

    // Read for the number of points
    char s[MAXLINE];
    if (fgets(s, MAXLINE, fp_mesh) == NULL)
    {
        fprintf(stderr, "Error reading file %s\n", MESH_FILE);
        exit(2);     
    }
    if (sscanf(s, "n_points %i", &n_points) != 1)
    {
        fprintf(stderr, "Did not find n_points in %s\n", MESH_FILE);
        exit(3);
    }

    // printf("n_points = %i\n", n_points);

    // Read for the points
    int i_file;
    mesh = malloc(sizeof(float) * n_points);
    if (mesh == NULL)
    {
        fprintf(stderr, "Error allocating memory for mesh\n");
        exit(5);
    }

    for (int i = 0; i < n_points; i++)
    {
        if (fgets(s, MAXLINE, fp_mesh) == NULL)
        {
            fprintf(stderr, "Error reading file %s. Missing points\n", MESH_FILE);
            exit(3);
        }

        if (sscanf(s, "%i %f", &i_file, &mesh[i]) != 2)
        {
            fprintf(stderr, "Error reading file %s. Missing points\n", MESH_FILE);
            exit(3);
        }

        if (i != i_file)
        {
            fprintf(stderr, "Error reading file %s. Wrong point sequence\n", MESH_FILE);
            exit(3);
        }
    }

    // Some cleanup
    fgets(s, MAXLINE, fp_mesh);
    if (!feof(fp_mesh))
    {
        fprintf(stderr, "Error reading file %s. EOF was not reached\n", MESH_FILE);
        exit(3);
    }
    fclose(fp_mesh);

    // printf("Done reading mesh!\n");

    // while (--n_points >= 0)
    // while (n_points-- > 0)
    /*for (int i = 0; i < n_points; i++)
    {
        printf("    %f\n", mesh[i]);
    }*/

    // print_matrix("Mesh:\n", mesh, n_points, 1);

   return;
}

// Reads settings from the SETTINGS_FILE
void read_settings(void)
{
    printf("Reading settings file ...\n");

    FILE *fp_settings = fopen(SETTINGS_FILE, "r");
    if (fp_settings == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", SETTINGS_FILE);
        exit(2);
    }

    int entries_read = fscanf(fp_settings, "Solver controls v%*i u_init %f phi_init %f phi_left %f phi_right %f rho %f gamma %f scheme %s",
    &settings.u_init, &settings.phi_init, &settings.phi_left, &settings.phi_right, &settings.rho, &settings.gamma, settings.scheme_name);

    fclose(fp_settings);

    if (entries_read != 7)
    {
        fprintf(stderr, "Did not find required entries in %s\n", SETTINGS_FILE);
        exit(3);
    }

    // Convert scheme name to scheme int
    if (!strcmp(settings.scheme_name, CDS_name))
    {
        settings.scheme_int = CDS;
    }
    else if (!strcmp(settings.scheme_name, UDS_name))
    {
        settings.scheme_int = UDS;
    }
    else
    {
        fprintf(stderr, "Wrong scheme!");
        exit(6);
    }

    // printf("Solver controls:\nu_init %f\nphi_init %f\nphi_left %f\nphi_right %f\nrho %f\ngamma %f\nscheme %s\n",
    // settings.u_init, settings.phi_init, settings.phi_left, settings.phi_right, settings.rho, settings.gamma, settings.scheme);

    return;
}

// Initialise the solution
void initialise(void)
{
    printf("Initializing the solution ...\n");

    u = malloc(sizeof(float) * n_points);
    if (u == NULL)
    {
        fprintf(stderr, "Error allocating memory for velocity\n");
        exit(5);
    }

    phi = malloc(sizeof(float) * n_points);
    if (phi == NULL)
    {
        fprintf(stderr, "Error allocating memory for phi\n");
        exit(5);
    }

    for (int i = 0; i < n_points; i++)
    {
        u[i] = settings.u_init;
        phi[i] = settings.phi_init;
    }

    phi[0] = settings.phi_left;
    phi[n_points-1] = settings.phi_right;

    // print_matrix("u:\n", u, n_points, 1);
    // print_matrix("phi:\n", phi, n_points, 1);

    n_pts_mtrx = n_points - 2;

    A = malloc(sizeof(float) * (int)pow(n_pts_mtrx, 2));
    if (A == NULL)
    {
        fprintf(stderr, "Error allocating memory for coefficient matrix\n");
        exit(5);
    }

    Q = malloc(sizeof(float) * (n_pts_mtrx));
    if (Q == NULL)
    {
        fprintf(stderr, "Error allocating memory for sources matrix\n");
        exit(5);
    }

    for (int i = 0; i < n_pts_mtrx; i++)
    {
        for (int j = 0; j < n_pts_mtrx; j++)
        {
            // A[j + i * (n_pts_mtrx)] = j + i * (n_pts_mtrx);
            A[j + i * n_pts_mtrx] = 0.0;
        }

        Q[i] = 0.0;
    }

    // print_matrix("A:\n", A, n_pts_mtrx, 2);
    // print_matrix("Q:\n", Q, n_pts_mtrx, 1);

    return;
}

// Construct matrices
void construct_matrices(void)
{
    printf("Constructing matrices ...\n");

    // Compute coefficients
    // The first node
    A[1] = A_e(1);
    A[0] = - (A_w(1) + A_e(1));
    Q[0] = - A_w(1) * phi[0];

    // The middle nodes
    for (int i = 2, i_mtrx = 0; i < n_pts_mtrx; i++)
    {
        // Find the corresponding index in the A matrix
        i_mtrx = (i - 1) * (n_pts_mtrx + 1);
        // A[i][i-1], i == i - 1
        A[i_mtrx - 1] = A_w(i);
        // A[i][i+1], i == i - 1
        A[i_mtrx + 1] = A_e(i);
        // A[i][i], i == i - 1
        A[i_mtrx] = - (A_w(i) + A_e(i));
        // Q[i], i == i - 1
        Q[i-1] = 0.0;
    }

    // The last node
    A[n_pts_mtrx * n_pts_mtrx - 2] = A_w(n_pts_mtrx);
    A[n_pts_mtrx * n_pts_mtrx - 1] = - (A_e(n_pts_mtrx) + A_w(n_pts_mtrx));
    Q[n_pts_mtrx - 1] = - A_e(n_pts_mtrx) * phi[n_pts_mtrx + 1];

    // print_matrix("A:\n", A, n_pts_mtrx, 2);
    // print_matrix("Q:\n", Q, n_pts_mtrx, 1);

    return;
}

// Solve the system using the TDMA
void solve_tdma(void)
{
    printf("Solving the system ...\n");

    // Forward elimination
    // [Re]calculate Ap and Q*
    float Ap[n_pts_mtrx];
    float Q_star[n_pts_mtrx];

    Ap[0] = A[0];
    Q_star[0] = Q[0];

    for (int i = 1, i_mtrx = 0; i < n_pts_mtrx; i++)
    {
        // Find the corresponding index in the A matrix
        i_mtrx = i * (n_pts_mtrx + 1);
        Ap[i] = A[i_mtrx] - A[i_mtrx - 1] * A[i_mtrx - n_pts_mtrx] / Ap[i - 1];
        Q_star[i] = Q[i] - A[i_mtrx - 1] * Q_star[i - 1] / Ap[i - 1];
    }

    // print_matrix("Ap:\n", Ap, n_pts_mtrx, 1);
    // print_matrix("Q_star:\n", Q_star, n_pts_mtrx, 1);

    // Back substitution (the phi matrix is n_points long)
    phi[n_pts_mtrx - 1 + 1] = Q_star[n_pts_mtrx - 1] / Ap[n_pts_mtrx - 1];

    for (int i = n_pts_mtrx - 2, i_mtrx = 0; i >= 0; i--)
    {
        i_mtrx = i * (n_pts_mtrx + 1);
        phi[i + 1] = (Q_star[i] - A[i_mtrx + 1] * phi[i + 1 + 1]) / Ap[i];
    }

    // print_matrix("phi:\n", phi, n_points, 1);

    return;
}

// Write results
void write_results(void)
{
    printf("Writing results ...\n");

    FILE *fp_results = fopen(RESULTS_FILE, "w");
    if (fp_results == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", RESULTS_FILE);
        exit(2);
    }

    fprint_matrix(fp_results, "phi:\n", phi, n_points, 1);

    fclose(fp_results);

    return;
}

// Cleanup
void free_mem(void)
{
    printf("Freeing memory ...\n");

    free(mesh);
    free(u);
    free(phi);
    free(A);
    free(Q);

    return;
}

// Print the supplied matrix beautifully
void print_matrix(char *s, float *M, int dim_size, int n_dims)
{
    printf("%s", s);

    int i = 0;
    int l = pow(dim_size, n_dims);

    while (i < l)
    {
        printf("%f", M[i]);

        if ((++i % dim_size) && (n_dims != 1))
        {
            printf("   ");
        }
        else
        {
            printf("\n");
        }
    }

    return;
}

// Print the supplied matrix beautifully to a file
void fprint_matrix(FILE *out, char *s, float *M, int dim_size, int n_dims)
{
    fprintf(out, "%s", s);

    int i = 0;
    int l = pow(dim_size, n_dims);

    while (i < l)
    {
        fprintf(out, "%f", M[i]);

        if ((++i % dim_size) && (n_dims != 1))
        {
            fprintf(out, "   ");
        }
        else
        {
            fprintf(out, "\n");
        }
    }

    return;
}

// Return the "west" coefficient
float A_w(int i)
{
    float a_w;

    // convection
    switch (settings.scheme_int)
    {
        case CDS:
            a_w = -settings.rho * u[i] / (mesh[i+1] - mesh[i-1]);
            break;

        case UDS:
            a_w = min(-settings.rho * u[i], 0.0) / (mesh[i] - mesh[i-1]);
            break;

        default:
            fprintf(stderr, "Wrong scheme!");
            exit(6);
    }
    
    // diffusion
    a_w += -2.0 * settings.gamma / (mesh[i+1] - mesh[i-1]) / (mesh[i] - mesh[i-1]);

    return a_w;
}

// Return the "east" coefficient
float A_e(int i)
{
    float a_e;

    // convection
    switch (settings.scheme_int)
    {
        case CDS:
            a_e = settings.rho * u[i] / (mesh[i+1] - mesh[i-1]);
            break;

        case UDS:
            a_e = min(settings.rho * u[i], 0.0) / (mesh[i+1] - mesh[i]);
            break;

        default:
            fprintf(stderr, "Wrong scheme!");
            exit(6);
    }

    // diffusion
    a_e += -2.0 * settings.gamma / (mesh[i+1] - mesh[i-1]) / (mesh[i+1] - mesh[i]);

    return a_e;
}

// Returns the smaller of the given arguments
float min(float a, float b)
{
    return (a < b) ? a : b;
}
