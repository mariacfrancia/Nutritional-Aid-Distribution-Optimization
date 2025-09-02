sets
    i poblados /Chipane,Saute,Camo-Camo,Faimba,Mazive,Macuene,Fernando,Banguene,Funhalouro/
    obj /Z2, Z3, Z4/;
    set pagos /distancia, coste, equidad/;
;

alias(i,j);
alias(pagos, p);

parameters
    d(i) demanda de la aldea i /Chipane 30,Saute 10,Camo-Camo 40,Faimba 40,Mazive 30,Macuene 20,Fernando 10,Banguene 40, Funhalouro 80/
    CS coste semanal de tener un centro abierto /50/
    CVV Coste de viaje de un vehículo en vacío /8/
    CTC Coste de transporte de carga por Kg /0.3/
    M cota superior /300/
    w1 peso distancia solución de compromiso apartado 4  /0.4/
    w2 peso coste solución de compromiso apartado 4 /0.2/
    w3 peso equidad solución de compromiso apartado 4 /0.4/
    p1 peso distancia solución de compromiso apartado 5 /0.393/
    p2 peso coste solución de compromiso apartado 5 /0.071/
    p3 peso equidad solución de compromiso apartado 5 /0.536/
    Z11 óptimo modelo 2 con 80%
    Z22 óptimo modelo 3 con 80%
    Z33 óptimo modelo 4 con 80%
    max_col(pagos) valor máximo de cada columna en matriz de pagos
    min_col(pagos) valor mínimo de cada columna en matriz de pagos
    paso2 paso de incremento para Z2
    paso3 paso de incremento para Z3
    paso4 paso de incremento para Z4

    eps2 cota superior para Z2
    eps3 cota superior para Z3
    eps4 cota superior para Z4
;

table
    t(i,j) tiempo que tardo en llegar de i a j
                 Chipane Saute Camo-Camo Faimba Mazive Macuene Fernando Banguene Funhalouro
    Chipane        0      4       3       99      99     99       3       99        3
    Saute          4      0       2       99      99     99       2       99       99
    Camo-Camo      3      2       0       99      99      2       3       99       99
    Faimba         99    99      99        0       3      3       99      99       99
    Mazive         99    99      99        3       0     99       99       2        3
    Macuene        99    99       2        3      99      0       99      99       99
    Fernando       3      2       3       99      99     99       0       99       99
    Banguene       99    99      99       99       2     99       99       0        2
    Funhalouro     3     99      99       99       3     99       99       2        0
;


variables
    X(i) si construimos centro en i
    Y(i,j) si i abastece a j
    R(i,j) cantidad que abastece i a j
    V(i,j) si trasporto carga de i a j
    Q(i,j) cantidad que transporto de i a j
    mdesv máxima desviación entre la demanda y la carga que llega
    MAXIMO para el caso de Tchebychev
    
    Z1 objetivo 1
    Z2 objetivo 2
    Z3 objetivo 3
    Z4 objetivo 4
    Z4_compromiso1 objetivo solución compromiso apartado 4
    Z4_compromisoinf objetivo solución compromiso apartado 4 con distancia Tchebychev
    Z5_compromiso objetivo solución compromiso apartado 5
;

binary variables X, Y, V;
positive variables R, Q, mdesv, cost, dist;

equations
    OBJETIVO1 objetivo 1
    OBJETIVO2 objetivo 2
    OBJETIVO22 objetivo 2 con el 80% de la demanda
    OBJETIVO3 objetivo 3
    OBJETIVO4 objetivo 4
    OBJETIVO4_compromiso1 objetivo 4 compromiso
    OBJETIVO4_compromiso_INF objetivo 4 compromiso distancia Tchebychev
    OBJETIVO5_compromiso objetivo 5 compromiso
    
    CAMINO_ABAST solo abastecemos por caminos posibles 
    ABAST_CENTRO solo abastecen centros
    RECUBRIMIENTO todos los polados deben ser abastecidos por un centro
    AUTOABAST todos los centros se autoabastecen
    MAX_CENTROS máximo 4 centros
    CAMINO_CARGA solo transportamos por caminos usados
    CAMINO_EXIST solo podemos transportar por caminos que existan
    NO_FUN a Funhalouro no entra carga
    FUN_SALIDA de Funhalouro sale toda la carga
    FLUJO de cada poblado lo que entra menos lo que sale es las demandas que abastece
    AUTOCAMINO no usamos caminos de i a i
    FUN_SALIDA80 de Funhalouro sale el 80% de la carga
    FLUJO80 de cada poblado lo que entra menos lo que sale es la cantidad real que abastece
    ABAST_POSIBLE solo abastecemos carga en conexiones posibles
    MAX_ABAST abastecemos como máximo la demanda
    MAX_DESV máxima desviación
    
    UNA_ENTRADA a cada poblado solo llega el camión una vez (evitamos ciclos de flujo)
    
    CAMINO_CARGA80 solo transportamos por caminos usados teniendo en cuenta 80% de carga
    TOTAL_R total de carga cuando es de un 80%
    RECUBRIMIENTO1 cada poblado debe ser abastecido 1 vez
    
    RESTR_EPS2 restriccion epsilon para Z2
    RESTR_EPS3 restriccion epsilon para Z3
    RESTR_EPS4 restriccion epsilon para Z4
    
    MAXIMA_DISTANCIA para el caso de Tchebychev en el modelo de compromiso4
    MAXIMO_COSTE para Tchebychev en el modelo de compromiso4
    MAXIMA_EQUIDAD para Tchebychev en el modelo de compromiso4
;
    
*funciones objetivo:

OBJETIVO1.. Z1 =E= sum(i, X(i));

OBJETIVO2.. Z2 =E= sum((i,j), 2*Y(i,j)*d(j)*t(i,j));

OBJETIVO22.. Z2 =E= sum((i,j), R(i,j)*t(i,j))*2;

OBJETIVO3.. Z3 =E= sum(i, X(i)*CS) + CVV*sum((i,j), V(i,j)*t(i,j))  + CTC*sum((i,j), Q(i,j)*t(i,j));

OBJETIVO4.. Z4 =E= mdesv;

OBJETIVO4_compromiso1.. Z4_compromiso1 =E= w1*(Z2-min_col("distancia"))/(max_col("distancia")-min_col("distancia")) +w2*(Z3-min_col("coste"))/(max_col("coste")-min_col("coste"))+w3*(Z4-min_col("equidad"))/(max_col("equidad")-min_col("equidad"));

OBJETIVO4_compromiso_INF.. Z4_compromisoinf =E= MAXIMO;

OBJETIVO5_compromiso.. Z5_compromiso =E= p1*(Z2-min_col("distancia"))/(max_col("distancia")-min_col("distancia")) +p2*(Z3-min_col("coste"))/(max_col("coste")-min_col("coste"))+p3*(Z4-min_col("equidad"))/(max_col("equidad")-min_col("equidad"));


*restricciones:

CAMINO_ABAST(i,j).. t(i,j)*Y(i,j) =L= 3*X(i);

ABAST_CENTRO(i,j).. Y(i,j) =L= X(i);

RECUBRIMIENTO1(j).. sum(i, Y(i,j)) =E= 1;
RECUBRIMIENTO(j).. sum(i, Y(i,j)) =L= 1;

AUTOABAST(i).. Y(i,i) =E= X(i);

MAX_CENTROS.. sum(i, X(i)) =L= 4;

CAMINO_CARGA(i,j).. Q(i,j) =L= V(i,j)*M;

CAMINO_CARGA80(i,j).. Q(i,j) =L= V(i,j)*M*0.8;

CAMINO_EXIST(i,j).. t(i,j)*V(i,j) =L= 5;

NO_FUN.. sum(i, V(i, 'Funhalouro')) =E= 0;

FUN_SALIDA.. sum(j, Q('Funhalouro',j)) =E= sum(i, d(i)) - sum(j, d(j)*Y('Funhalouro',j));

FLUJO(j)$(not sameas(j,'Funhalouro')).. sum(i, Q(i,j)) - sum(i, Q(j,i)) =E= sum(i, d(i)*Y(j,i));

AUTOCAMINO(i).. V(i,i) =E= 0;

FUN_SALIDA80.. sum(j, Q('Funhalouro',j)) =E= 0.8*sum(i, d(i)) - sum(j, R('Funhalouro',j));

FLUJO80(j)$(not sameas(j,'Funhalouro')).. sum(i, Q(i,j)) - sum(i, Q(j,i)) =E= sum(i, R(j,i));

ABAST_POSIBLE(i,j).. R(i,j) =L= Y(i,j)*M;

MAX_ABAST(j).. sum(i, R(i,j)) =L= d(j);

MAX_DESV(i).. (d(i) - sum(j, R(j,i))) / d(i) =L= mdesv;

UNA_ENTRADA(j).. sum(i, V(i,j)) =L= 1;

TOTAL_R.. sum((i,j),R(i,j)) =E= M*0.8;

* ε-Restricciones: minimizando la funcion objetivo 
RESTR_EPS2.. Z2 =L= eps2;
RESTR_EPS3.. Z3 =L= eps3;
RESTR_EPS4.. Z4 =L= eps4;

MAXIMA_DISTANCIA.. MAXIMO =G= w1*(Z2-min_col("distancia"))/(max_col("distancia")-min_col("distancia"));

MAXIMO_COSTE.. MAXIMO =G= w2*(Z3-min_col("coste"))/(max_col("coste")-min_col("coste"));

MAXIMA_EQUIDAD.. MAXIMO =G= w3*(Z4-min_col("equidad"))/(max_col("equidad")-min_col("equidad"));



*====================================================================

*modelo de optimizacion
MODEL MODELO1 /OBJETIVO1, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO1, AUTOABAST/;

*MIP porque es un problema mixto con variables continuas y binarias
SOLVE MODELO1 MINIMIZING Z1 USING MIP;
DISPLAY X.l, Y.l, Z1.l;

*====================================================================
*Modelo de distancia:
MODEL MODELO2 /OBJETIVO2, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO1, AUTOABAST, MAX_CENTROS/

*MIP porque es un problema mixto con variables continuas y binarias
SOLVE MODELO2 MINIMIZING Z2 USING MIP;
DISPLAY X.l, Y.l, Z2.l;

*====================================================================
*Modelo de coste:
MODEL MODELO3 /OBJETIVO3, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO1, AUTOABAST, MAX_CENTROS, CAMINO_CARGA, CAMINO_EXIST, NO_FUN, FUN_SALIDA
                FLUJO, AUTOCAMINO, UNA_ENTRADA/

*MIP porque es un problema mixto con variables continuas y binarias
SOLVE MODELO3 MINIMIZING Z3 USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z3.l;

*====================================================================

*Modelo de equidad:
MODEL MODELO4 /OBJETIVO4, MAX_ABAST, MAX_DESV,TOTAL_R/
                
*MIP porque es un problema mixto con variables continuas y binarias
SOLVE MODELO4 MINIMIZING Z4 USING MIP;
DISPLAY Z4.l, R.l;


*====================================================================
*                      MATRIZ DE PAGOS                              *
*====================================================================  
* ε-Restricciones.
*hemos planteado 2 formas de hacer la matriz de pagos:
parameter matriz_pagos(pagos,pagos);
parameter matriz_pagos2(pagos,pagos);

*-----------------------------------------------------------------------------------
*fijamos Z2 (distancia): -----------------------------------------------------------
MODEL MODELO411 /OBJETIVO22, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R/
SOLVE MODELO411 MINIMIZING Z2 USING MIP;
DISPLAY Z2.l, X.l, Y.l, Q.l, V.l, R.l;
parameter Z11;
Z11 = Z2.l;
matriz_pagos('distancia','distancia') = Z11;
matriz_pagos2('distancia','distancia') = Z11;
equation OBJETIVO1_1;
OBJETIVO1_1.. Z11 =E= sum((i,j), 2*R(i,j)*t(i,j));

* hallamos Z3 en funcion de ese Z2: (MATRIZ DE PAGOS 2)
MODEL MODELO412 /OBJETIVO3, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R,OBJETIVO1_1/
SOLVE MODELO412 MINIMIZING Z3 USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z3.l, R.l;
parameter Z12;
Z12 = Z3.l;
matriz_pagos2('distancia','coste') = Z12;
equation OBJETIVO1_2;
OBJETIVO1_2.. Z12 =E= sum(i, X(i)*CS) + CVV*sum((i,j), V(i,j)*t(i,j))  + CTC*sum((i,j), Q(i,j)*t(i,j));

** hallamos Z4 en funcion de ese Z2: (MATRIZ DE PAGOS 2)
MODEL MODELO4133 /OBJETIVO4,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R,OBJETIVO1_1/
SOLVE MODELO4133 MINIMIZING Z4 USING MIP;
DISPLAY Z4.l,R.l;
matriz_pagos2('distancia','equidad') = Z4.l;

** hallamos Z4 en funcion de ese Z2 y ese Z1: (MATRIZ DE PAGOS)
MODEL MODELO413 /OBJETIVO4,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R,OBJETIVO1_1,OBJETIVO1_2/
SOLVE MODELO413 MINIMIZING Z4 USING MIP;
DISPLAY Z4.l,R.l;
matriz_pagos('distancia','equidad') = Z4.l;
parameter Z13;
Z13 = Z4.l;
equation OBJETIVO1_3;
OBJETIVO1_3(i).. Z13 =G= (d(i) - sum(j, R(j,i))) / d(i);

*hallamos Z3 en funcion de ese Z2 y este ultimo Z4: (MATRIZ DE PAGOS)
MODEL MODELO4122 /OBJETIVO3, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R,OBJETIVO1_1,OBJETIVO1_3/
SOLVE MODELO4122 MINIMIZING Z3 USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z3.l, R.l;
matriz_pagos('distancia','coste') = Z3.l;

*-------------------------------------------------------------------------------
*fijamos Z3 (coste): -----------------------------------------------------------
MODEL MODELO422 /OBJETIVO3, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R/
SOLVE MODELO422 MINIMIZING Z3 USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z3.l, R.l;
Z22 = Z3.l;
matriz_pagos('coste','coste') = Z22;
matriz_pagos2('coste','coste') = Z22;
equation OBJETIVO2_2;
OBJETIVO2_2.. Z22 =E= sum(i, X(i)*CS) + CVV*sum((i,j), V(i,j)*t(i,j))  + CTC*sum((i,j), Q(i,j)*t(i,j));

* hallamos Z2 en funcion de ese Z3: (MATRIZ DE PAGOS 2)
MODEL MODELO421 /OBJETIVO22, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO2_2/
SOLVE MODELO421 MINIMIZING Z2 USING MIP;
DISPLAY Z2.l, X.l, Y.l, Q.l, V.l, R.l;
matriz_pagos2('coste','distancia') = Z2.l;
parameter Z21;
Z21 = Z2.l;
equation OBJETIVO2_1;
OBJETIVO2_1.. Z21 =E= sum((i,j), 2*R(i,j)*t(i,j));

** hallamos Z4 en funcion de ese Z3: (MATRIZ DE PAGOS 2)
MODEL MODELO4233 /OBJETIVO4,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO2_2/
SOLVE MODELO4233 MINIMIZING Z4 USING MIP;
DISPLAY Z4.l,R.l;
matriz_pagos2('coste','equidad') = Z4.l;

** hallamos Z2 en funcion de ese Z3 y ese Z2: (MATRIZ DE PAGOS)
MODEL MODELO423 /OBJETIVO4,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO2_2, OBJETIVO2_1/
SOLVE MODELO423 MINIMIZING Z4 USING MIP;
DISPLAY Z4.l,R.l;
matriz_pagos('coste','equidad') = Z4.l;
parameter Z23;
Z23 = Z4.l;
equation OBJETIVO2_3;
OBJETIVO2_3(i).. Z23 =G= (d(i) - sum(j, R(j,i))) / d(i);

* hallamos Z2 en funcion de ese Z3 y este ultimo Z4: (MATRIZ DE PAGOS)
MODEL MODELO4211 /OBJETIVO22, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO2_2, OBJETIVO2_3/
SOLVE MODELO4211 MINIMIZING Z2 USING MIP;
DISPLAY Z2.l, X.l, Y.l, Q.l, V.l, R.l;
matriz_pagos('coste','distancia') = Z2.l;

*---------------------------------------------------------------------------------
*fijamos Z4 (equidad): -----------------------------------------------------------
MODEL MODELO433 /OBJETIVO4,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R/
SOLVE MODELO433 MINIMIZING Z4 USING MIP;
DISPLAY Z4.l,R.l;
Z33 = Z4.l;
matriz_pagos('equidad','equidad') = Z33;
matriz_pagos2('equidad','equidad') = Z33;
equation OBJETIVO3_3;
OBJETIVO3_3(i).. Z33 =G= (d(i) - sum(j, R(j,i))) / d(i);

* hallamos Z2 en funcion de ese Z4: (MATRIZ DE PAGOS 2)
MODEL MODELO431 /OBJETIVO22,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO3_3/
SOLVE MODELO431 MINIMIZING Z2 USING MIP;
DISPLAY Z2.l, X.l, Y.l, Q.l, V.l, R.l;
matriz_pagos2('equidad','distancia') = Z2.l;
parameter Z31;
Z31 = Z2.l;
equation OBJETIVO3_1;
OBJETIVO3_1.. Z31 =E= sum((i,j), 2*R(i,j)*t(i,j));

** hallamos Z3 en funcion de ese Z4: (MATRIZ DE PAGOS 2)
MODEL MODELO4322 /OBJETIVO3, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO3_3/
SOLVE MODELO4322 MINIMIZING Z3 USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z3.l, R.l;
matriz_pagos2('equidad','coste') = Z3.l;

** hallamos Z3 en funcion de ese Z4: (MATRIZ DE PAGOS)
MODEL MODELO432 /OBJETIVO3, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO3_3, OBJETIVO3_1/
SOLVE MODELO432 MINIMIZING Z3 USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z3.l, R.l;
matriz_pagos('equidad','coste') = Z3.l;
parameter Z32;
Z32 = Z3.l;
equation OBJETIVO3_2;
OBJETIVO3_2.. Z32 =E= sum(i, X(i)*CS) + CVV*sum((i,j), V(i,j)*t(i,j))  + CTC*sum((i,j), Q(i,j)*t(i,j));

* hallamos Z2 en funcion de ese Z4 y este ultimo Z2: (MATRIZ DE PAGOS)
MODEL MODELO4311 /OBJETIVO22,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO3_3, OBJETIVO3_2/
SOLVE MODELO4311 MINIMIZING Z2 USING MIP;
DISPLAY Z2.l, X.l, Y.l, Q.l, V.l, R.l;
matriz_pagos('equidad','distancia') = Z2.l;


display matriz_pagos;
display matriz_pagos2;


*====================================================================
*                   FRONTERAS DE PARETO                             *
*====================================================================
*lo hemos planteado de 2 maneras:

loop(pagos,
  max_col(pagos) = smax(p, matriz_pagos2(p,pagos));
  min_col(pagos) = smin(p, matriz_pagos2(p,pagos));
);

set g iteraciones /0*10/;
alias(g,kg)

set
  pz /Z2, Z3, Z4/;

parameter resultados(g, g, pz);

paso2 =  (max_col("distancia")-min_col("distancia"))/10;
paso3 =  (max_col("coste")-min_col("coste"))/10;
paso4 =  (max_col("equidad")-min_col("equidad"))/10;


* 1- minimizando 1 criterio con 1 restriccion: ------------------------------------------------------------------------------------
* es decir, 2 a 2:

** Para distancia y coste: función objetivo distancia y el objetivo restringido coste ==========================================
MODEL modeloPareto2 /CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO22, OBJETIVO3,RESTR_EPS3/;
parameter resultados2(g, pz);
parameter epsilon3(g);
loop(g, 
    eps3 = min_col("coste") + paso3 * (ord(g) - 1);

    SOLVE modeloPareto2 MINIMIZING Z2 USING MIP;

    resultados2(g, 'Z2') = Z2.l;
    resultados2(g, 'Z3') = Z3.l;
    
    epsilon3(g) = eps3
    
);


** Para distancia y coste: función objetivo coste y el objetivo restringido distancia
MODEL modeloPareto3 /CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO22, OBJETIVO3,RESTR_EPS2/;
parameter resultados3(g, pz);
loop(g, 
    eps2 = min_col("distancia") + paso2 * (ord(g) - 1);

    SOLVE modeloPareto3 MINIMIZING Z3 USING MIP;

    resultados3(g, 'Z2') = Z2.l;
    resultados3(g, 'Z3') = Z3.l;
);

** Para distancia y equidad: función objetivo distancia y el objetivo restringido equidad
MODEL modeloPareto4 /CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO22, OBJETIVO4,RESTR_EPS4/;
parameter resultados4(g, pz);
loop(g, 
    eps4 = min_col("equidad") + paso4 * (ord(g) - 1);

    SOLVE modeloPareto4 MINIMIZING Z2 USING MIP;

    resultados4(g, 'Z2') = Z2.l;
    resultados4(g, 'Z4') = Z4.l;
);

** Para distancia y equidad: función objetivo equidad y el objetivo restringido distancia
MODEL modeloPareto5 /CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO22, OBJETIVO4,RESTR_EPS2/;
parameter resultados5(g, pz);
loop(g, 
    eps2 = min_col("distancia") + paso2 * (ord(g) - 1);

    SOLVE modeloPareto5 MINIMIZING Z4 USING MIP;

    resultados5(g, 'Z2') = Z2.l;
    resultados5(g, 'Z4') = Z4.l;
);

** Para coste y equidad: función objetivo coste y el objetivo restringido equidad
MODEL modeloPareto6 /CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO3, OBJETIVO4,RESTR_EPS4/;
parameter resultados6(g, pz);
loop(g, 
    eps4 = min_col("equidad") + paso4 * (ord(g) - 1);

    SOLVE modeloPareto6 MINIMIZING Z3 USING MIP;

    resultados6(g, 'Z3') = Z3.l;
    resultados6(g, 'Z4') = Z4.l;
);

** Para coste y equidad: función objetivo equidad y el objetivo restringido coste
MODEL modeloPareto7 /CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO3, OBJETIVO4,RESTR_EPS3/;
parameter resultados7(g, pz);
loop(g, 
    eps3 = min_col("coste") + paso3 * (ord(g) - 1);

    SOLVE modeloPareto7 MINIMIZING Z4 USING MIP;

    resultados7(g, 'Z3') = Z3.l;
    resultados7(g, 'Z4') = Z4.l;
);


DISPLAY resultados2, resultados3, resultados4, resultados5, resultados6,resultados7

* 2- minimizando 1 criterio con 2 restricciones: ------------------------------------------------------------------------------------
*usando los 3 criterios a la vez:

* Resolver modelo minimizando Z4, con Z2 y Z3 como restricciones
MODEL modeloPareto /OBJETIVO4,CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R, OBJETIVO22, OBJETIVO3,RESTR_EPS2,RESTR_EPS3/;
 
loop(g,
  loop(kg,
    eps2 = min_col("distancia") + paso2 * (ord(g) - 1);
    eps3 = min_col("coste")     + paso3 * (ord(kg) - 1);

    SOLVE modeloPareto MINIMIZING Z4 USING MIP;

    resultados(g, kg, 'Z4') = Z4.l;
    resultados(g, kg, 'Z2') = Z2.l;
    resultados(g, kg, 'Z3') = Z3.l;
  );
);

DISPLAY resultados


*====================================================================
*                 MODELOS DE COMPROMISO                             *
*====================================================================

MODEL MODELO4_COMPROMISO_1 /OBJETIVO4_compromiso1, OBJETIVO22, OBJETIVO3,OBJETIVO4, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R/
SOLVE MODELO4_COMPROMISO_1 MINIMIZING Z4_compromiso1 USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z4_compromiso1.l, R.l, Z2.l,Z3.l,Z4.l;

MODEL MODELO4_COMPROMISO_INF /OBJETIVO4_compromiso_INF,OBJETIVO22, OBJETIVO3,OBJETIVO4, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R,MAXIMA_DISTANCIA,MAXIMO_COSTE,MAXIMA_EQUIDAD/
SOLVE MODELO4_COMPROMISO_INF MINIMIZING Z4_compromisoinf USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z4_compromisoinf.l, R.l, Z2.l,Z3.l,Z4.l;

*====================================================================

MODEL MODELO5_COMPROMISO /OBJETIVO5_compromiso, OBJETIVO22, OBJETIVO3,OBJETIVO4, CAMINO_ABAST, ABAST_CENTRO, RECUBRIMIENTO, AUTOABAST, MAX_CENTROS, CAMINO_CARGA80, CAMINO_EXIST, NO_FUN
                FUN_SALIDA80, FLUJO80, AUTOCAMINO, ABAST_POSIBLE, MAX_ABAST, MAX_DESV, UNA_ENTRADA,TOTAL_R/

SOLVE MODELO5_COMPROMISO MINIMIZING Z5_compromiso USING MIP;
DISPLAY X.l, Y.l, Q.l, V.l, Z5_compromiso.l, R.l, Z2.l,Z3.l,Z4.l;