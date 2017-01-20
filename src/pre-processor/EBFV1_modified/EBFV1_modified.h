/*
 * EBFV1_modified.h
 *
 *  Created on: Oct 17, 2014
 *      Author: rogerio
 */

#ifndef EBFV1_MODIFIED_H_
#define EBFV1_MODIFIED_H_

#include "EBFV1__pre-processors.h"
#include "GeomData.h"
#include "rockProp.h"

void EBFV1_modified_preprocessor_2D(pMesh theMesh, GeomData* pGCData);
void assemblyMatrix_A(const double*, const double*, const double*, const double*, const double*, const double*,const double*, const double*, const double*,double*,double*, double*);
void calculateGeomCoefficients(pEntity, double *, double *, double &);
void calculateMatrix_F(pEntity, const double*, const double*, double, double*, double*, double*);
void calculeMatrix_E(pEntity face, GeomData* pGCData, const double* K, const double* Cij, double *E_ij);
void calculateMatrix_E(pEntity, GeomData*, const double*, const double*, double*, double*, double*);
void calculeMatrix_G(pEntity edge, GeomData* pGCData, const double* K, const double* Cij, double *G_ij);
void calculateMatrix_G(pEntity, GeomData*, const double*, const double*, double*, double*, double*);

#endif /* EBFV1_MODIFIED_H_ */

/*------------------------------------------------------------------------------------------------------------
Outubro/2014

Deducao da equacao discreta do termo de pressao para uma formulacao de volumes finitos baseadas em arestas


Malhas 2-D de triangulos:
-----------------------------------------------------------------------------------------
Gij__2x2 = -(K*lamb*Cij*Lij/|Lij|)*[1 -1; -1 1];

Eij__2x4=-(K*lamb/2)*[I-Lij*Lij' I-Lij*Lij'; Lij*Lij'-I Lij*Lij'-I];

Fij_omega__4x2 = Cij/2*[1/Vi 1/Vi; -1/Vj -1/Vj];

Fij_gamma__4x2 = [5/6Dij/Vi Dij/Vi; Dij/Vj 5/6Dij/Vj];



Equacao de pressão (difusão/eliptica)
-div(K*lamb*grad_P) = q

Integrando ambos os lados em omega:
E*F*u+G*u = Q



Dado uma malha formada por um unico elemento triangular (arestas: I-J, I-K, J-K)
-----------------------------------------------------------------------------------------
I-J
Gij = -(K*lamb*Cij*Lij/|Lij|)*[1 -1; -1 1];
Eij=-(K*lamb/2)*[I-Lij*Lij' I-Lij*Lij'; Lij*Lij'-I Lij*Lij'-I];
Fij_omega = Cij/2*[1/Vi 1/Vi; -1/Vj -1/Vj];
Fij_gamma = [5/6Dij/Vi Dij/Vi; Dij/Vj 5/6Dij/Vj];

I-K
Gik = -(K*lamb*Cik*Lik/|Lik|)*[1 -1; -1 1];
Eik=-(K*lamb/2)*[I-Lik*Lik' I-Lik*Lik'; Lik*Lik'-I Lik*Lik'-I];
Fik_omega = Cik/2*[1/Vi 1/Vi; -1/Vk -1/Vk];
Fik_gamma = [5/6Dik/Vi Dik/Vi; Dik/Vk 5/6Dik/Vk];

J-K
Gjk = -(K*lamb*Cjk*Ljk/|Ljk|)*[1 -1; -1 1];
Ejk=-(K*lamb/2)*[I-Ljk*Ljk' I-Ljk*Ljk'; Ljk*Ljk'-I Ljk*Ljk'-I];
Fjk_omega = Cjk/2*[1/Vj 1/Vj; -1/Vk -1/Vk];
Fjk_gamma = [5/6Djk/Vj Djk/Vj; Djk/Vk 5/6Djk/Vk];

Montando a matrix E global: E_npxnp.dim => E_3x6

                        I                                         J                            K
    |E_ij[0,0]+E_ik[0,0] E_ij[0,1]+E_ik[0,1]          E_ij[0,2]         E_ij[0,3]                E_ik[0,2]           E_ik[0,3]     |
E = |     E_ij[1,0]           E_ij[1,1]         E_ij[1,2]+E_jk[0,0] E_ij[1,3]+E_jk[0,1]           E_jk[0,2]           E_jk[0,3]     |
    |     E_ik[1,0]           E_ik[1,1]               E_jk[1,0]         E_jk[1,1]           E_ik[1,2]+E_jk[1,2] E_ik[1,3]*E_jk[1,3]|


Montando a matrix F global: F_np.dimxnp => F_6x3

    |F_ij[0,0]+F_ik[0,0]      F_ij[0,1]           F_ik[0,1]     |
    |F_ij[1,0]+F_ik[1,0]      F_ij[1,1]           F_ik[1,1]     |
F = |     F_ij[2,0]      F_ij[2,1]+F_jk[0,0]      F_jk[0,1]     |
    |     F_ij[3,0]      F_ij[3,1]+F_jk[1,0]      F_jk[1,1]     |
    |     F_ik[2,0]           F_jk[2,0]      F_ik[2,1]+F_jk[2,1]|
    |     F_ik[3,0]           F_jk[3,0]      F_ik[3,1]+F_jk[3,1]|

Montando a matriz G global: G_npxnp => G_3x3

    |G_ij[0,0]+G_ik[0,0]      G_ij[0,1]           G_ik[0,1]     |
G = |     G_ij[1,0]      G_ij[1,1]+G_jk[0,0]      G_jk[0,1]     |
    |     G_ik[1,0]           G_jk[1,0]      G_ik[1,1]+G_jk[1,1]|

Montando a matriz global do elemento: A = E*F+G
       |A11 A12 A13|
A_ge = |A21 A22 A23|
       |A31 A32 A33|

    A11 = (E_ij[0,0]+E_ik[0,0])*(F_ij[0,0]+F_ik[0,0])+
          (E_ij[0,1]+E_ik[0,1])*(F_ij[1,0]+F_ik[1,0])+
           E_ij[0,2]*F_ij[2,0]+
           E_ij[0,3]*F_ij[3,0]+
           E_ik[0,2]*F_ik[2,0]+
           E_ik[0,3]*F_ik[3,0]+
           G_ij[0,0]+G_ik[0,0]


    A12 = (E_ij[0,0]+E_ik[0,0])*F_ij[0,1]+
          (E_ij[0,1]+E_ik[0,1])*F_ij[1,1]+
           E_ij[0,2]*(F_ij[2,1]+F_jk[0,0])+
           E_ij[0,3]*(F_ij[3,1]+F_jk[1,0])+
           E_ik[0,2]*F_jk[2,0]+
           E_ik[0,3]*F_jk[3,0]+
           G_ij[0,1]

    A13 = (E_ij[0,0]+E_ik[0,0])*F_ik[0,1]+
          (E_ij[0,1]+E_ik[0,1])*F_ik[1,1]+
           E_ij[0,2]*F_jk[0,1])+
           E_ij[0,3]*F_jk[1,1])+
           E_ik[0,2]*(F_ik[2,1]+F_jk[2,1])+
           E_ik[0,3]*(F_ik[3,1]+F_jk[3,1])+
           G_ik[0,1]

    A21 =  E_ij[1,0]*(F_ij[0,0]+F_ik[0,0])+
           E_ij[1,1]*(F_ij[1,0]+F_ik[1,0])+
           (E_ij[1,2]+E_jk[0,0])*F_ij[2,0]+
           (E_ij[1,3]+E_jk[0,1])*F_ij[3,0]+
           E_jk[0,2]*F_ik[2,0]+
           E_jk[0,3]*F_ik[3,0]+
           G_ij[1,0]

    A22 =  E_ij[1,0]*F_ij[0,1]+
           E_ij[1,1]*F_ij[1,1]+
           (E_ij[1,2]+E_jk[0,0])*(F_ij[2,1]+F_jk[0,0])+
           (E_ij[1,3]+E_jk[0,1])*(F_ij[3,1]+F_jk[1,0])+
           E_jk[0,2]*F_jk[2,0]+
           E_jk[0,3]*F_jk[3,0]+
           G_ij[1,1]+G_jk[0,0]

    A23 =  E_ij[1,0]*F_ik[0,1]+
           E_ij[1,1]*F_ik[1,1]+
           (E_ij[1,2]+E_jk[0,0])*F_jk[0,1]+
           (E_ij[1,3]+E_jk[0,1])*F_jk[1,1]+
           E_jk[0,2]*(F_ik[2,1]+F_jk[2,1])+
           E_jk[0,3]*(F_ik[3,1]+F_jk[3,1])+
           G_jk[0,1]

    A31 = E_ik[1,0]*(F_ij[0,0]+F_ik[0,0])+
          E_ik[1,1]*(F_ij[1,0]+F_ik[1,0])+
          E_jk[1,0]*F_ij[2,0]+
          E_jk[1,1]*F_ij[3,0]+
          (E_ik[1,2]+E_jk[1,2])*F_ik[2,0]+
          (E_ik[1,3]*E_jk[1,3])*F_ik[3,0]+
          G_ik[1,0]

    A32 = E_ik[1,0]*F_ij[0,1]+
          E_ik[1,1]*F_ij[1,1]+
          E_jk[1,0]*(F_ij[2,1]+F_jk[0,0])+
          E_jk[1,1]*(F_ij[3,1]+F_jk[1,0])+
          (E_ik[1,2]+E_jk[1,2])*F_jk[2,0]+
          (E_ik[1,3]*E_jk[1,3])*F_jk[3,0]+
          G_jk[1,0]

    A33 = E_ik[1,0]*F_ik[0,1]+
          E_ik[1,1]*F_ik[1,1]+
          E_jk[1,0]*F_jk[0,1]+
          E_jk[1,1]*F_jk[1,1]+
          (E_ik[1,2]+E_jk[1,2])*(F_ik[2,1]+F_jk[2,1])+
          (E_ik[1,3]*E_jk[1,3])*(F_ik[3,1]+F_jk[3,1])+
          G_ik[1,1]+G_jk[1,1]

Desejamos agora escrever uma matriz por aresta. Esta por sua vez será uma matriz 2x3. Assim, a matriz global final do problema será
montada por aresta. Mas antes de ser adicionada a global final, a matriz por aresta precisa ser multiplicada pela mobilidade calculada
sobre a aresta e isto depende apenas das saturações dos nós da aresta.

    A_ij[0,0] = E_ij[0,0]*(F_ij[0,0]+F_ik[0,0])+E_ij[0,1]*(F_ij[1,0]+F_ik[1,0])+E_ij[0,2]*F_ij[2,0]+E_ij[0,3]*F_ij[3,0]+G_ij[0,0];
    A_ij[0,1] = E_ij[0,0]*F_ij[0,1]+E_ij[0,1]*F_ij[1,1]+E_ij[0,2]*(F_ij[2,1]+F_jk[0,0])+E_ij[0,3]*(F_ij[3,1]+F_jk[1,0])+G_ij[0,1];
    A_ij[0,2] = E_ij[0,0]*F_ik[0,1]+E_ij[0,1]*F_ik[1,1]+E_ij[0,2]*F_jk[0,1]+E_ij[0,3]*F_jk[1,1];
    A_ij[1,0] = E_ij[1,0]*(F_ij[0,0]+F_ik[0,0])+E_ij[1,1]*(F_ij[1,0]+F_ik[1,0])+E_ij[1,2]*F_ij[2,0]+E_ij[1,3]*F_ij[3,0]+G_ij[1,0];
    A_ij[1,1] = E_ij[1,0]*F_ij[0,1]+E_ij[1,1]*F_ij[1,1]+E_ij[1,2]*(F_ij[2,1]+F_jk[0,0])+E_ij[1,3]*(F_ij[3,1]+F_jk[1,0])+G_ij[1,1];
    A_ij[1,2] = E_ij[1,0]*F_ik[0,1]+E_ij[1,1]*F_ik[1,1]+E_ij[1,2]*F_jk[0,1]+E_ij[1,3]*F_jk[1,1];

    A_jk[0,0] = E_jk[0,0]*F_ij[2,0]+E_jk[0,1]*F_ij[3,0]+E_jk[0,2]*F_ik[2,0]+E_jk[0,3]*F_ik[3,0];
    A_jk[0,1] = E_jk[0,0]*(F_ij[2,1]+F_jk[0,0])+E_jk[0,1]*(F_ij[3,1]+F_jk[1,0])+E_jk[0,2]*F_jk[2,0]+E_jk[0,3]*F_jk[3,0]+G_jk[0,0];
    A_jk[0,2] = E_jk[0,0]*F_jk[0,1]+E_jk[0,1]*F_jk[1,1]+E_jk[0,2]*(F_ik[2,1]+F_jk[2,1])+E_jk[0,3]*(F_ik[3,1]+F_jk[3,1])+G_jk[0,1];
    A_jk[1,0] = E_jk[1,0]*F_ij[2,0]+E_jk[1,1]*F_ij[3,0]+E_jk[1,2]*F_ik[2,0]+E_jk[1,3]*F_ik[3,0];
    A_jk[1,1] = E_jk[1,0]*(F_ij[2,1]+F_jk[0,0])+E_jk[1,1]*(F_ij[3,1]+F_jk[1,0])+E_jk[1,2]*F_jk[2,0]+E_jk[1,3]*F_jk[3,0]+G_jk[1,0];
    A_jk[1,2] = E_jk[1,0]*F_jk[0,1]+E_jk[1,1]*F_jk[1,1]+E_jk[1,2]*(F_ik[2,1]+F_jk[2,1])+E_jk[1,3]*(F_ik[3,1]+F_jk[3,1])+G_jk[1,1];

    A_ik[0,0] = E_ik[0,0]*(F_ij[0,0]+F_ik[0,0])+E_ik[0,1]*(F_ij[1,0]+F_ik[1,0])+E_ik[0,2]*F_ik[2,0]+E_ik[0,3]*F_ik[3,0]+G_ik[0,0];
    A_ik[0,1] = E_ik[0,0]*F_ij[0,1]+E_ik[0,1])*F_ij[1,1]+E_ik[0,2]*F_jk[2,0]+E_ik[0,3]*F_jk[3,0];
    A_ik[0,2] = E_ik[0,0]*F_ik[0,1]+E_ik[0,1]*F_ik[1,1]+E_ik[0,2]*(F_ik[2,1]+F_jk[2,1])+E_ik[0,3]*(F_ik[3,1]+F_jk[3,1])+G_ik[0,1];
    A_ik[1,0] = E_ik[1,0]*(F_ij[0,0]+F_ik[0,0])+E_ik[1,1]*(F_ij[1,0]+F_ik[1,0])+E_ik[1,2]*F_ik[2,0]+E_ik[1,3]*F_ik[3,0]+G_ik[1,0];
    A_ik[1,1] = E_ik[1,0]*F_ij[0,1]+E_ik[1,1]*F_ij[1,1]+E_ik[1,2]*F_jk[2,0]+E_ik[1,3]*F_jk[3,0];
    A_ik[1,2] = E_ik[1,0]*F_ik[0,1]+E_ik[1,1]*F_ik[1,1]+E_ik[1,2]*(F_ik[2,1]+F_jk[2,1])+E_ik[1,3]*(F_ik[3,1]+F_jk[3,1])+G_ik[1,1];

As matrizes acima (A_ij, A_jk, A_ik) são matrizes por aresta de um mesmo elemento. Como uma aresta pode ser compartilhada por dois elementos,
a cada uma delas será adicionada a matriz correspondente do elemento vizinho. Isso se dá da seguinte forma:

A_ij = A_ij_E1 + A_ij_E2;

onde A_ij_E1 e A_ij_E2 são matrizes de elementos adjacente, ou seja, que compartilham a memsma aresta IJ. O mesmo raciocínio se dá para
as demais matrizes.

Todas as matrizes são funções de E e G. A contribuicao de F é vista como uma constante, pois nao depende de Sw.

Na formulacao original proposto de volumes finitos, a matriz global para um problema com multidominios é dada por:

A = (E*F+G)_dom1 + (E*F+G)_dom2 + ... + (E*F+G)_domN

Se imaginarmos cada elemento da malha como um dominio, podemos escrever a equacao acima como:

A = A_ge_dom1 + A_ge_dom2 + ... + A_ge_domN

Cada matriz global do elemento que representa uma dominio K (A_ge_domk) pode ser decomposta em matrizes de arestas do elemento, entao basta
um loop sobre as arestas da malha para que a matriz global final seja montada.




Exibindo EBFV1_modificada.txt…

 */
