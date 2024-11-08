#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 512

float V1[N], V2[N], V3[N], V4[N];
float Mat[N][N], MatDD[N][N];

void InitData(){
    int i,j;
    srand(334411);
    for( i = 0; i < N; i++ )
        for( j = 0; j < N; j++ ){
            Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
            if ( (abs(i - j) <= 3) && (i != j))
            MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
            else if ( i == j )
            MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
            else MatDD[i][j] = 0.0;
        }

        for( i = 0; i < N; i++ ){
            V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
            V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
            V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
        }
}

void PrintVect( float vect[N], int from, int numel ){
    for(  int i = from; i < from + numel; i++){
        printf("%f ", vect[i]);
    }
}

void PrintRow( float mat[N][N], int row, int from, int numel ){
    for(int i = from; i < from + numel; i++){
        printf("%f ", mat[row][i]);  
    }
}

void MultEscalar( float vect[N], float vectres[N], float alfa ){
    for(int i = 0; i < N; i++){
        vectres[i] = vect[i] * alfa;
    }
}

float Scalar (float vect1[N], float vect2[N]){
    float resultado = 0;
    for(int i = 0; i < N; i++){
        resultado +=  vect1[i] * vect2[i];
    }
    return resultado;
}

float Magnitude( float vect[N] ){
    float sum = 0;
    for (int i = 0; i < N; i++) {
        sum += (vect[i] * vect[i]);
    }
    return sqrt(sum);
}

int Ortogonal( float vect1[N], float vect2[N]){
    float resultado = Scalar(vect1, vect2);
    if(resultado == 0){
        return 1; //Si que son ortogonals
    } else {
        return 0; //NO son ortogonals
    }
}

void Projection( float vect1[N], float vect2[N], float vectres[N] ){
    float producte_escalar = Scalar(vect1, vect2);
    float magnitud_v = Magnitude(vect2);

    for (int i = 0; i < N; i++) {
        vectres[i] = (producte_escalar / magnitud_v) * vect2[i];
    }
}

float Infininorm( float M[N][N] ){
    float maxim = 0;
    for(int i = 0; i < N; i++){
        float suma = 0;
        for(int j = 0; j < N; j++){
            suma = suma + fabs(M[i][j]);
        }
        if(suma > maxim){
            maxim = suma;
        }      
    }
    return maxim;
}

float Onenorm( float M[N][N] ){
    float maxim = 0;
    for(int i = 0; i < N; i++){
        float suma = 0;
        for(int j = 0; j < N; j++){
            suma = suma + fabs(M[j][i]);
        }
        if(suma > maxim){
            maxim = suma;
        }      
    }
    return maxim;
}

float NormFrobenius( float M[N][N] ){
    float suma = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            suma = suma + (M[i][j] * M[i][j]);
        }
    }
    float resultado = sqrt(suma);
    return resultado;
}

int DiagonalDom( float M[N][N] ){
    float diagonal;
    float suma;
    for (int i = 0; i < N; i++) {
        suma = 0;
        diagonal = fabs(M[i][i]);
        for(int j = 0; j < N; j++){
            if (i != j){
                suma = suma + fabs(M[i][j]);
            }
        }
        if (diagonal < suma) {
            return 0; // No es diagonal dominante
        }
    }
    return 1; // Si es diagonal dominante
}

void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ){
    for(int i = 0; i < N; i++){
        vectres[i] = 0;
        for(int j = 0; j < N; j++){
            vectres[i] = vectres[i] + (M[i][j] * vect[j]);
        }
    }
}

int Jacobi(float A[N][N], float b[N], float x_result[N], unsigned max_iter) {
    int aplicable = 1;
    float x_antiguo[N];
    float x_nuevo[N];
    float suma;

    for (int i = 0; i < N; i++) {
        x_antiguo[i] = 0;
        x_nuevo[i] = 0;
    }
    
    if (DiagonalDom(A) == 1) {
        for (int paso = 0; paso < max_iter; paso++) {
            for (int i = 0; i < N; i++) {
                suma = b[i];
                for (int j = 0; j < N; j++) {
                    if (j != i) {
                        suma -= A[i][j] * x_antiguo[j];
                    }
                }
                x_nuevo[i] = suma / A[i][i];
            }
            for (int i = 0; i < N; i++) {
                x_antiguo[i] = x_nuevo[i];
            }
        }
        for (int i = 0; i < N; i++) {
            x_result[i] = x_antiguo[i];
        }
    } else {
        aplicable = 0;
    }

    return aplicable;
}


int main(){
    //Inicializamos vectores y matrius
    InitData();

    //A
    printf("V1 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V1, 0, 10);
    PrintVect(V1, 256, 10);
    printf("\nV2 del 0 al 9 i del 256 al 265:\n");
    PrintVect(V2, 0, 10);
    PrintVect(V2, 256, 10);
    printf("\nV3 del 0 al 9 i del 256 al 265::\n");
    PrintVect(V3, 0, 10);
    PrintVect(V3, 256, 10);

    //B
    printf("\n\nMat fila 0 i fila 100 del 0 al 9:\n");
    PrintRow(Mat, 0, 0, 10);
    PrintRow(Mat, 100, 0, 10);

    //C
    printf("\n\nMatDD fila 0 del 0 al 9 i fila 100 del 95 al 104:\n");
    PrintRow(MatDD, 0, 0, 10);
    PrintRow(MatDD, 100, 95, 10);

    //D
    //MAT
    float infininormMAT = Infininorm(Mat);
    printf("\n\nInfininorma de Mat = %f\n", infininormMAT);
    float onenormMAT = Onenorm(Mat);
    printf("\nNorma ú de Mat = %f\n", onenormMAT);
    float frobeniusMAT = NormFrobenius(Mat);
    printf("\nNorma de Frobenius de Mat = %f\n", frobeniusMAT);
    if (DiagonalDom(Mat) == 1){
        printf("\nLa matriu Mat és diagonal dominant\n");
    } else {
        printf("\nLa matriu Mat no és diagonal dominant\n");

    }

    //MATDD
    float infininormMATDD = Infininorm(MatDD);
    printf("\n\nInfininorma de MatDD = %f\n", infininormMATDD);
    float onenormMATDD = Onenorm(MatDD);
    printf("\nNorma ú de MatDD = %f\n", onenormMATDD);
    float frobeniusMATDD = NormFrobenius(MatDD);
    printf("\nNorma de Frobenius de MatDD = %f\n", frobeniusMATDD);
    if (DiagonalDom(MatDD) == 1){
        printf("\nLa matriu MatDD és diagonal dominant\n");
    } else {
        printf("\nLa matriu MatDD no és diagonal dominant\n");

    }

    //E
    //A V1 * V2
    float producteV1V2 = Scalar(V1, V2);
    printf("\n\nEscalar <V1,V2> = %f", producteV1V2);
    //A V1 * V3
    float producteV1V3 = Scalar(V1, V3);
    printf("Escalar <V1,V3> = %f", producteV1V3);
    //A V2 * V3
    float producteV2V3 = Scalar(V2, V3);
    printf("Escalar <V2,V3> =  %f", producteV2V3);

    //F
    float magnitudV1 = Magnitude(V1);
    float magnitudV2 = Magnitude(V2);
    float magnitudV3 = Magnitude(V3);
    printf("\n\nMagnitud V1,V2 i V3 = %f %f %f", magnitudV1, magnitudV2, magnitudV3);

    //G
    if (Ortogonal(V1, V2) == 1){
        printf("\n\nV1 i V2 són ortogonals\n");
    }
    if (Ortogonal(V1, V3) == 1){
        printf("\nV1 i V3 són ortogonals\n");
    }
    if (Ortogonal(V2, V3) == 1){
        printf("\nV2 i V3 són ortogonals\n");
    }

    //H
    MultEscalar(V3, V4, 2);
    printf("\n\nEls elements 0 al 9 i 256 al 265 del resultat de multiplicar V3x2.0 són: \n");
    PrintVect(V4, 0, 10);
    PrintVect(V4, 256, 10);

    //I 
    //ERROR AQUI
    Projection(V2, V3, V4);
    printf("\n\nEls elements 0 a 9 del resultat de la projecció de V2 sobre V3 són: \n");
    PrintVect(V4, 0, 10);

    Projection(V1, V2, V4);
    printf("\n\nEls elements 0 a 9 del resultat de la projecció de V1 sobre V2 són: \n");
    PrintVect(V4, 0, 10);

    //J
    Matriu_x_Vector(Mat, V2, V4);
    printf("\n\nEls elements 0 a 9 del resultat de la multiplicació de Mat per v2 són: \n");
    PrintVect(V4, 0, 10);

    //K
    //MatDD*X = V3
    // 1 iter MatDD*X = V3
    int aplicableMatDD_1 = Jacobi(MatDD, V3, V4, 1);
    if (aplicableMatDD_1 == 1) {
        printf("\nEls elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
        PrintVect(V4, 0, 10);
    }

    // 1000 iter MatDD*X = V3
    int aplicableMatDD_1000 = Jacobi(MatDD, V3, V4, 1000);
    if (aplicableMatDD_1000 == 1) {
        printf("\nEls elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
        PrintVect(V4, 0, 10);
    }

    if (aplicableMatDD_1 == 0 || aplicableMatDD_1000 == 0) {
        printf("\nLa matriu MDD no és diagonal dominant, no es pot aplicar Jacobi\n");
    }

    //Mat*X = V3
    // 1 iter Mat*X = V3
    int aplicableMat_1 = Jacobi(Mat, V3, V4, 1);
    if (aplicableMat_1 == 1) {
        printf("\nEls elements 0 a 9 de la solució (1 iter)  del sistema d'equacions són:\n");
        PrintVect(V4, 0, 10);
    } 

    // 1000 iter Mat*X = V3
    int aplicableMat_1000 = Jacobi(Mat, V3, V4, 1000);
    if (aplicableMat_1000 == 1) {
        printf("\nEls elements 0 a 9 de la solució (1000 iters)  del sistema d'equacions són:\n");
        PrintVect(V4, 0, 10);
    }

    if (aplicableMat_1 == 0 || aplicableMat_1000 == 0) {
        printf("\nLa matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
    }

}
