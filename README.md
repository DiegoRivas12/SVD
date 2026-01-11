# Descomposición en Valores Singulares (SVD) con LAPACK

Este repositorio contiene un ejemplo de cálculo de la Descomposición en Valores Singulares (SVD) de una matriz 2×2 utilizando LAPACK en C, demostrando que un ejemplo teórico proporcionado contiene errores en sus cálculos.

---

## Contenido del Repositorio

- `svd_example.c` — Programa en C que calcula la SVD con LAPACK
- `ejecucion_lapack.jpg` — Captura de pantalla de la ejecución
- `README.md` — Este documento

---

## Tabla de Contenidos

1. [Compilación y Ejecución](#compilación-y-ejecución)
2. [Resultado de la Ejecución](#resultado-de-la-ejecución)
3. [Cálculo Teórico Correcto](#cálculo-teórico-correcto)
4. [Comparación: Valores Teóricos vs LAPACK](#comparación-valores-teóricos-vs-lapack)
5. [Conclusiones](#conclusiones)

---

## Compilación y Ejecución

### Requisitos
- Compilador GCC
- Biblioteca LAPACK instalada
- Biblioteca BLAS instalada

### Instalación de dependencias

**Ubuntu/Debian:**
```bash
sudo apt-get install liblapack-dev libblas-dev
```

**macOS:**
```bash
brew install lapack openblas
```

### Compilar
```bash
gcc -o svd_example svd_example.c -llapack -lblas -lm
```

### Ejecutar
```bash
./svd_example
```

---

## Resultado de la Ejecución

![Ejecución LAPACK](ejecucion_lapack.jpg)

### Interpretación de los Resultados

El programa calcula la SVD de la matriz:
```
A = [ 1.0   -0.8 ]
    [ 0.0    1.0 ]
```

**Salida del programa:**

1. **Valores Singulares (Σ)**
```
   σ₁ = 1.477033
   σ₂ = 0.677033
```

2. **Matriz U** (vectores singulares izquierdos)
```
   [ -0.828067   0.560629 ]
   [  0.560629   0.828067 ]
```

3. **Matriz V^T** (vectores singulares derechos transpuestos)
```
   [ -0.560629   0.828067 ]
   [  0.828067   0.560629 ]
```

4. **Verificación**: A reconstruida = U × Σ × V^T
```
   [ 1.000000  -0.800000 ]
   [ 0.000000   1.000000 ]
```

✓ La reconstrucción es exacta (error < 10⁻¹⁵)

---

## Cálculo Teórico Correcto

### Objetivo
Demostrar que los valores singulares del ejemplo proporcionado (σ₁ = 1.62, σ₂ = 1.00) son **incorrectos**, y calcular los valores correctos que coinciden con LAPACK.

---

### Paso 1: Cálculo de A^T × A
```
A^T = [ 1.0    0.0 ]
      [-0.8    1.0 ]

A^T × A = [ 1.00   -0.80 ]
          [-0.80    1.64 ]
```

---

### Paso 2: Ecuación Característica
```
det(A^T × A - λI) = 0

(1 - λ)(1.64 - λ) - 0.64 = 0
λ² - 2.64λ + 1.00 = 0
```

---

### Paso 3: Autovalores (Fórmula Cuadrática)
```
λ = (2.64 ± √(2.64² - 4)) / 2
λ = (2.64 ± √2.9696) / 2
λ = (2.64 ± 1.7233) / 2

λ₁ = 2.1817
λ₂ = 0.4583
```

---

### Paso 4: Valores Singulares
```
σ₁ = √2.1817 = 1.4771 ✓
σ₂ = √0.4583 = 0.6770 ✓
```

**Comparación:**
```
Ejemplo dado (INCORRECTO):  σ₁ = 1.62,   σ₂ = 1.00
Cálculo correcto:           σ₁ = 1.4771, σ₂ = 0.6770
LAPACK:                     σ₁ = 1.4770, σ₂ = 0.6770
```

---

### Paso 5: Demostración del Error

---

### Paso 5.1: Generación de las Matrices del Ejemplo Incorrecto

Para entender de dónde provienen las matrices del ejemplo dado, debemos **forzar** los valores singulares incorrectos y calcular las matrices U y V correspondientes.

#### Valores singulares forzados del ejemplo:
```
σ₁ = 1.62  (requiere λ₁ = 2.62)
σ₂ = 1.00  (requiere λ₂ = 1.00)
```

#### Cálculo de vectores con valores incorrectos

**Para λ₁ = 2.62 (FORZADO):**

Resolvemos (A^T × A - 2.62·I)v₁ = 0:
```
[ 1.00-2.62   -0.80     ] [x]   [0]
[-0.80        1.64-2.62 ] [y] = [0]

[-1.62   -0.80] [x]   [0]
[-0.80   -0.98] [y] = [0]
```

De la segunda ecuación:
```
-0.80x - 0.98y = 0
x = -0.98y / 0.80 = -1.225y
```

Tomando y = 1:
```
v₁_no_normalizado = [-1.225]
                    [ 1.000]

||v₁|| = √(1.225² + 1²) = √2.5006 = 1.581

v₁ = [-1.225/1.581]   [-0.78]
     [ 1.000/1.581] = [ 0.62]  ← Valores del ejemplo
```

**Para λ₂ = 1.00 (FORZADO):**

Resolvemos (A^T × A - 1.00·I)v₂ = 0:
```
[ 0.00   -0.80] [x]   [0]
[-0.80    0.64] [y] = [0]
```

De la primera ecuación:
```
-0.80y = 0  →  pero esto da y = 0 (inconsistente)
```

De la segunda ecuación:
```
-0.80x + 0.64y = 0
x = 0.80y
```

Tomando y = -1 (ajuste de signo para coincidir con ejemplo):
```
v₂_no_normalizado = [-0.80]
                    [-1.00]

||v₂|| = √(0.64 + 1) = √1.64 = 1.281

v₂ = [-0.80/1.281]   [-0.62]
     [-1.00/1.281] = [-0.78]  ← Valores del ejemplo
```

#### Matriz V del ejemplo (con valores forzados):
```
V = [-0.78   -0.62]
    [ 0.62   -0.78]

V^T = [-0.78    0.62]
      [-0.62   -0.78]  ← Coincide con el ejemplo dado
```

#### Cálculo de U con valores singulares incorrectos

Usando uᵢ = (1/σᵢ) × A × vᵢ con los σ incorrectos:

**Para u₁:**
```
A × v₁ = [ 1.0   -0.8] × [-0.78]
         [ 0.0    1.0]   [ 0.62]

       = [-1.276]
         [ 0.62]

u₁ = (1/1.62) × [-1.276]   [-0.79]
                [ 0.62 ] = [ 0.38]  ← Valores del ejemplo
```

**Para u₂:**
```
A × v₂ = [ 1.0   -0.8] × [-0.62]
         [ 0.0    1.0]   [-0.78]

       = [ 0.00]
         [-0.78]

u₂ = (1/1.00) × [ 0.00]   [ 0.00]
                [-0.78] = [-0.78]  ← Valores del ejemplo
```

#### Matrices del ejemplo incorrecto:
```
U = [-0.79    0.00]
    [ 0.38   -0.78]  ← Coincide con el ejemplo dado

Σ = [1.62   0.00]
    [0.00   1.00]

V^T = [-0.78    0.62]
      [-0.62   -0.78]
```

#### Verificación con valores incorrectos:

Si multiplicamos estas matrices:
```
U × Σ × V^T = [-0.79    0.00] × [1.62   0.00] × [-0.78    0.62]
              [ 0.38   -0.78]   [0.00   1.00]   [-0.62   -0.78]
```

**El resultado NO reconstruye exactamente la matriz A original**, porque los valores singulares no son los correctos matemáticamente.

---

### Resumen del error en cascada:
```
Autovalores incorrectos
        ↓
Valores singulares incorrectos  
        ↓
Vectores V incorrectos
        ↓
Vectores U incorrectos
        ↓
Reconstrucción A imprecisa
```

**Conclusión**: Las matrices del ejemplo (U, Σ, V^T) solo existen porque se **forzaron** valores singulares que no son solución de la ecuación característica de A^T·A. Aunque las matrices tienen la estructura correcta de una SVD, no representan la verdadera descomposición de la matriz A dada.

---

**Verificación de los valores incorrectos:**
```
Para λ = 2.62:
(2.62)² - 2.64(2.62) + 1.00 = 0.9476 ≠ 0 ✗

Para λ = 1.00:
(1.00)² - 2.64(1.00) + 1.00 = -0.64 ≠ 0 ✗
```

**Verificación de los valores correctos:**
```
Para λ = 2.1817:
(2.1817)² - 2.64(2.1817) + 1.00 ≈ 0 ✓

Para λ = 0.4583:
(0.4583)² - 2.64(0.4583) + 1.00 ≈ 0 ✓
```

**Los valores del ejemplo NO satisfacen la ecuación característica.**

---

### Paso 6: Vectores Singulares Derechos (V)

**Para λ₁ = 2.1817:**
```
(A^T·A - λ₁I)v₁ = 0
-1.1817x - 0.80y = 0  →  x = -0.6770y

v₁ normalizado = [-0.5606]
                 [ 0.8281]
```

**Para λ₂ = 0.4583:**
```
(A^T·A - λ₂I)v₂ = 0
0.5417x - 0.80y = 0  →  x = 1.4771y

v₂ normalizado = [ 0.8281]
                 [ 0.5606]
```
```
V^T = [-0.5606   0.8281]
      [ 0.8281   0.5606]
```

---

### Paso 7: Vectores Singulares Izquierdos (U)

Usando uᵢ = (1/σᵢ) × A × vᵢ:
```
u₁ = [-0.8281]    u₂ = [ 0.5606]
     [ 0.5606]         [ 0.8281]

U = [-0.8281   0.5606]
    [ 0.5606   0.8281]
```

---

### Paso 8: Verificación
```
A = U × Σ × V^T = [ 1.0000  -0.8000] ✓
                  [ 0.0000   1.0000]
```

---

## Comparación: Valores Teóricos vs LAPACK

### Valores Singulares

| Valor singular | Ejemplo dado | Correcto (LAPACK) | Diferencia |
|---------------|--------------|-------------------|------------|
| σ₁            | 1.6200       | 1.4770            | 0.1430     |
| σ₂            | 1.0000       | 0.6770            | 0.3230     |

### Matriz U

| Componente | Ejemplo dado | LAPACK   |
|------------|--------------|----------|
| u₁₁        | -0.7900      | -0.8281  |
| u₂₁        | -0.3800      |  0.5606  |
| u₁₂        |  0.0000      |  0.5606  |
| u₂₂        | -0.7800      |  0.8281  |

### Matriz V^T

| Componente | Ejemplo dado | LAPACK   |
|------------|--------------|----------|
| v₁₁        | -0.7800      | -0.5606  |
| v₁₂        |  0.6200      |  0.8281  |
| v₂₁        | -0.6200      |  0.8281  |
| v₂₂        | -0.7800      |  0.5606  |

---

## Conclusiones

### 1. El ejemplo teórico contiene errores fundamentales

- **Autovalores incorrectos**: λ₁ = 2.62 y λ₂ = 1.00 no satisfacen la ecuación característica, pero eran necesarios para poder reconstruir las matrices que pedia el ejemplo.
- **Valores singulares incorrectos**: σ₁ = 1.62 y σ₂ = 1.00 son consecuencia del error anterior
- **Vectores incorrectos**: Todos los vectores U y V derivados son inválidos

### 2. LAPACK proporciona los valores correctos

- Valores singulares: σ₁ = 1.4770, σ₂ = 0.6770
- Verificados mediante cálculo analítico independiente
- Precisión numérica del orden de 10⁻¹⁵

### 3. Método de verificación

Este trabajo demuestra la importancia de:
- Verificar resultados teóricos con implementaciones numéricas estables
- Comprobar que los autovalores satisfacen la ecuación característica
- Usar bibliotecas especializadas (LAPACK) para cálculos críticos

### 4. Notas adicionales

- Los signos de vectores singulares tienen ambigüedad (v y -v son válidos)
- LAPACK usa algoritmos optimizados para estabilidad numérica
- La reconstrucción A = U·Σ·V^T es la prueba definitiva de corrección

---

## Referencias

- LAPACK Users' Guide: https://www.netlib.org/lapack/
- Documentación de LAPACK dgesvd: https://www.netlib.org/lapack/explore-html/
