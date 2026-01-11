#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Declaración de dgesvd de LAPACK
extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
                    double *a, int *lda, double *s,
                    double *u, int *ldu, double *vt, int *ldvt,
                    double *work, int *lwork, int *info);

int main() {
    int m = 2, n = 2;
    int lda = m, ldu = m, ldvt = n;

    // Matriz original A en column-major
    double A[4] = {1.0, 0.0,   // columna 1
                   -0.8, 1.0}; // columna 2
    double Acopy[4] = {1.0, 0.0, -0.8, 1.0}; // copia para LAPACK

    double U[4], VT[4], S[2];
    int info;
    int lwork = 5*(m>n? m:n);
    double *work = (double*)malloc(lwork*sizeof(double));

    char jobu = 'A', jobvt = 'A';

    // Llamada a LAPACK
    dgesvd_(&jobu, &jobvt, &m, &n, Acopy, &lda, S, U, &ldu, VT, &ldvt, work, &lwork, &info);

    if(info != 0) {
        printf("Error en dgesvd_, info = %d\n", info);
        free(work);
        return -1;
    }

    // --- Ajuste opcional de signos para coincidir con teoría ---
    // Multiplicamos columna 0 de U y fila 0 de VT por -1
    for(int i=0;i<m;i++) U[i + 0*ldu] *= -1;
    for(int j=0;j<n;j++) VT[0 + j*ldvt] *= -1;

    // --- Imprimir matriz original ---
    printf("Matriz original A:\n");
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++) printf("%f ", A[i + j*lda]);
        printf("\n");
    }
    printf("\n");

    // --- Imprimir valores singulares ---
    printf("Valores singulares (Sigma):\n");
    for(int i=0;i<n;i++) printf("%f ", S[i]);
    printf("\n\n");

    // --- Imprimir matriz U ---
    printf("Matriz U:\n");
    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++) printf("%f ", U[i + j*ldu]);
        printf("\n");
    }
    printf("\n");

    // --- Imprimir matriz V^T ---
    printf("Matriz V^T:\n");
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++) printf("%f ", VT[i + j*ldvt]);
        printf("\n");
    }
    printf("\n");

    // --- Multiplicación correcta paso a paso Arec = U*S*VT ---
    double Arec[4];
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            Arec[i + j*lda] = 0.0;
            for(int k=0;k<n;k++){  // n valores singulares
                Arec[i + j*lda] += U[i + k*ldu] * S[k] * VT[k + j*ldvt];
            }
        }
    }

    printf("Matriz A reconstruida (U*S*V^T):\n");
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++) printf("%f ", Arec[i + j*lda]);
        printf("\n");
    }
    printf("\n");

    free(work);
    return 0;
}
