# Test Suite para extract_genes_final.py

Este directorio contiene una suite completa de pruebas usando `pytest` para validar el funcionamiento del programa `extract_genes_final.py`.

## Estructura de Pruebas

Las pruebas están organizadas en las siguientes clases:

### 1. **TestReverseComplement** (6 pruebas)
Pruebas para la función `_reverse_complement()`:
- Secuencias básicas
- Nucleótidos ambiguos (N)
- Conversión de mayúsculas/minúsculas
- Secuencias largas
- Nucleótidos individuales

### 2. **TestValidatePath** (4 pruebas)
Pruebas para la validación de rutas de archivos:
- Archivos FASTA válidos
- Archivos GFF válidos
- Archivos inexistentes
- Extensiones incorrectas

### 3. **TestLoadFasta** (10 pruebas)
Pruebas para la carga de archivos FASTA:
- Secuencias individuales
- Múltiples secuencias
- Conversión a mayúsculas
- Secuencias multilínea
- Manejo de líneas vacías
- IDs duplicados
- Archivos vacíos (errores)
- Encabezados malformados (errores)
- Encabezados vacíos (errores)
- Encabezados con espacios

### 4. **TestParseGff** (11 pruebas)
Pruebas para el análisis de archivos GFF:
- Análisis básico
- Campos de características
- Coordenadas como enteros
- Hebra positiva/negativa
- Filtrado de características no-gen
- Análisis de atributos
- Fallback de gene_id
- Manejo de errores (columnas malformadas)
- Manejo de errores (coordenadas inválidas)
- Ignorar comentarios

### 5. **TestExtractGeneSeqs** (8 pruebas)
Pruebas para la extracción de secuencias de genes:
- Extracción básica
- Coordenadas correctas
- Hebra directa (+)
- Hebra inversa (-) con complemento inverso
- Errores de seqid inexistente
- Errores de coordenadas fuera de rango
- Manejo de genes sin nombre
- Información de hebra preservada

### 6. **TestWriteFastaOutput** (5 pruebas)
Pruebas para la escritura del archivo FASTA de salida:
- Escritura básica
- Formato de encabezado
- Múltiples genes
- Envoltura de líneas
- Creación de directorios

### 7. **TestIntegration** (2 pruebas)
Pruebas de integración del flujo completo:
- Flujo completo de extracción
- Flujo con genes en ambas hebras

## Ejecutar las Pruebas

### Ejecutar todas las pruebas:
```bash
pytest test_extract_genes/ -v
```

### Ejecutar pruebas de una clase específica:
```bash
pytest test_extract_genes/test_extract_genes.py::TestLoadFasta -v
```

### Ejecutar una prueba específica:
```bash
pytest test_extract_genes/test_extract_genes.py::TestLoadFasta::test_load_single_sequence -v
```

### Ejecutar con cobertura de código:
```bash
pytest test_extract_genes/ --cov=src --cov-report=html
```

### Ejecutar con salida detallada:
```bash
pytest test_extract_genes/ -vv -s
```

## Fixtures

El archivo `conftest.py` proporciona las siguientes fixtures reutilizables:

- `temp_dir`: Directorio temporal para archivos de prueba
- `sample_fasta`: Archivo FASTA de ejemplo válido
- `sample_gff`: Archivo GFF de ejemplo válido
- `malformed_fasta`: Archivo FASTA malformado
- `empty_header_fasta`: FASTA con encabezado vacío
- `duplicate_ids_fasta`: FASTA con IDs duplicados
- `malformed_gff`: Archivo GFF malformado
- `invalid_coordinates_gff`: GFF con coordenadas inválidas
- `nonexistent_seqid_gff`: GFF con seqid inexistente

## Resultados de Pruebas

Actualmente, todas las **46 pruebas pasan** con éxito:

```
=================================== 46 passed in 1.97s ====================================
```

## Notas

- Las pruebas usan fixtures de `pytest` para crear archivos temporales
- Se cubre tanto el comportamiento correcto como los casos de error
- Las pruebas de integración validan el flujo completo del programa
- Todos los archivos de prueba son autolimpiables (se usan directorios temporales)
