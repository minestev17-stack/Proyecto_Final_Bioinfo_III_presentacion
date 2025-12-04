'''
#Uso
python extract_genes.py --gff genes.gff --fasta genome.fasta --output genes.fna


Objetivo
Leer un archivo FASTA con el genoma completo.
Leer un archivo GFF.
Extraer las secuencias DNA correspondientes a las features gene.
Guardarlas en un archivo FASTA de salida.


Ejemplo de salida esperada (genes.fna)
>araC  gene_coords=3456-41020 strand=+
ATGCGTAGCTAGCTAGCTAGCTAA
>crp  gene_coords=3456-41020 strand=-
ATTTGCGCGGCGCGCGTTAG


Parte B — Extensión
Agregar:
--min-length N
Ejemplo:

python extract_genes.py --gff genes.gff --fasta genome.fasta --output genes.fna --min-length 300
Solo genes con longitud ≥ 300 serán exportados.

Requisitos técnicos
argparse obligatorio.
Funciones: load_fasta(), parse_gff(), extract_gene_seqs().
Manejo de errores con excepciones.
Docstrings y PEP8.
Pruebas con asserts + documento de pruebas, bien con pytest.
'''

import argparse
#TODO: importar otros módulos necesarios

def parse_arguments():
    '''
    TODO: Implementar el parseo de argumentos.
    --fasta: ruta al archivo FASTA del genoma. Es obligatorio. Tiene que terminar en .fasta o .fa o .fna
    --gff: ruta al archivo GFF. Es obligatorio. Tiene que terminar en .gff o .gff3
    --output: ruta al archivo FASTA de salida. Es obligatorio. Tiene que terminar en .fasta o .fa o .fna
    '''
    pass

def load_fasta(fasta_path): -> dict: [str, str]
    ''' 
    TODO: Implementar la carga del archivo FASTA. Validar formato.
    '''
    pass

def parse_gff(gff_path): -> list[dict]:
    '''
    TODO: Implementar el parseo del archivo GFF. Validar formato.
    '''
    pass

def write_fasta_output(output_path, gene_seqs_form_extracted):
    '''
    TODO: Implementar la escritura del archivo FASTA de salida.
    ''' 
    pass

def extract_gene_seqs(fasta_dict, gff_features): -> list[dict]:
    ''' 
    TODO: Implementar la extracción de secuencias de genes.
    '''
    
def main():
    #TODO: Ejecutar las funciones implementadas
    pass


if __name__ == '__main__':
    main()