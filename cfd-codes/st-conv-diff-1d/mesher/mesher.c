/*
 * Mesher for a steady convection-diffusion problem in 1D.
 * Version: 0.1
 * gcc -o mesher .\mesher.c
*/

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Requires
#define MESH_CONTROLS "meshControls.txt"
#define MESH_RESULTS "meshPoints.txt"
#define LEFT 0.0
#define RIGHT 1.0
#define N_POINTS 11
#define MAXLINE 20

// Prototypes
int readMeshControls(char *filename, float *left, float *right, int *n_points, float *ref_fac);
float sum_pow(float r, int n);

// Do stuff
int main(int argc, char *argv[])
{
    printf("Starting mesher ...\n");

    // Check inputs
    float left, right, ref_fac;
    int n_points;

    if (argc == 1)
    {
        printf("No file was given as an argument. Looking for %s ...\n", MESH_CONTROLS);
        readMeshControls(MESH_CONTROLS, &left, &right, &n_points, &ref_fac);
    }
    else if (argc == 2)
    {
        printf("Looking for mesh controls in %s ...\n", argv[1]);
        readMeshControls(argv[1], &left, &right, &n_points, &ref_fac);
    }
    else
    {
        fprintf(stderr, "Usage: .\\mesher.exe [meshControls]\n");
        exit(1);
    }
    
    // Constructing mesh
    printf("Constructing 1D mesh from x = %f to x = %f using %i points and refinement factor %f ...\n", left, right, n_points, ref_fac);

    // Prepare to write mesh
    printf("Writing mesh to %s\n", MESH_RESULTS);
    FILE *fp_mesh = fopen(MESH_RESULTS, "w");
    if (fp_mesh == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", MESH_RESULTS);
        exit(2);
    }
    fprintf(fp_mesh, "n_points %i\n", n_points);
    
    double dxmin = (right - left) / sum_pow(ref_fac, n_points - 1);
    // printf("    dxmin = %f\n", dxmin);

    double xi = left;
    for (int i = 0; i < n_points; i++)
    {
        xi = left + dxmin * sum_pow(ref_fac, i);
        // printf("    x%i = %f\n", i, xi);
        fprintf(fp_mesh, "%i %f\n", i, xi);
    }

    // Housekeeping
    if (ferror(fp_mesh))
    {
        fprintf(stderr, "Error writing mesh file\n");
        exit(4);
    }
    fclose(fp_mesh);

    printf("Done!\n");
    if (ferror(stdout))
    {
        fprintf(stderr, "Error writing stdout\n");
        exit(4);
    }

    return 0;
}

// Read for mesh controls in a supplied file
int readMeshControls(char *filename, float *left, float *right, int *n_points, float *ref_fac)
{
    printf("Reading %s ...\n", filename);
    
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(2);
    }

    char s[MAXLINE];
    
    fgets(s, MAXLINE, fp);
    if (sscanf(s, "left %f", left) != 1)
    {
        fprintf(stderr, "Missing the 'left' entry\n");
        exit(3);
    }
    
    fgets(s, MAXLINE, fp);
    if (sscanf(s, "right %f", right) != 1)
    {
        fprintf(stderr, "Missing the 'right' entry\n");
        exit(3);
    }
    
    fgets(s, MAXLINE, fp);
    if (sscanf(s, "n_points %d", n_points) != 1)
    {
        fprintf(stderr, "Missing the 'n_points' entry\n");
        exit(3);
    }
    
    fgets(s, MAXLINE, fp);
    if (sscanf(s, "ref_fac %f", ref_fac) != 1)
    {
        fprintf(stderr, "Missing the 'ref_fac' entry\n");
        exit(3);
    }

    fclose(fp);
    return 0;
}

// Return the sum of powers of the refinement factor r, i.e.
// sum from i = 0 to i = n of r**i
float sum_pow(float r, int n)
{
    // Base case
    if (n <= 0)
    {
        return 0;
    }
    
    // Recursion
    return pow(r, n) + sum_pow(r, n - 1);
}
