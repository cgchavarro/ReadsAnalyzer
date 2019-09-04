package uniandes.algorithms.readsanalyzer;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

/**
 * Simple script that simulates error free reads from a text in fasta format
 * @author Jorge Duitama
 *
 */
public class SimpleReadsSimulator {
	/**
	 * Main class that executes the program
	 * @param args Array of arguments:
	 * args[0]: Source sequence in fasta format. If many sequences are present, it only takes the first sequence
	 * args[1]: Length of the reads to simulate
	 * args[2]: Number of reads to simulate
	 * args[3]: Path to the output file
	 * args[4]: Number of errors in a args[5] number of bp
	 * @throws Exception If the fasta file can not be loaded
	 */
	public static void main(String[] args) throws Exception {
		String filename = args[0];
		int readLength = Integer.parseInt(args[1]);
		int numReads = Integer.parseInt(args[2]);
		String outFile = args[3];
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.setSequenceType(StringBuilder.class);
		QualifiedSequenceList sequences = handler.loadSequences(filename);
		if(sequences.size()==0) throw new Exception("No sequences found in file: "+filename);
		QualifiedSequence seq = sequences.get(0);
		String sequence = seq.getCharacters().toString();
		int seqLength = sequence.length();
		System.out.println("Length of the sequence to simulate reads: "+seqLength);
		double averageRD = ((double)numReads*readLength)/seqLength;
		System.out.println("Expected average RD: "+averageRD);
		char [] fixedQS = new char [readLength];
		Arrays.fill(fixedQS, '5');
		String fixedQSStr = new String(fixedQS);
		Random random = new Random();


		try (PrintStream out = new PrintStream(outFile)){
			// TODO: Generar lecturas aleatorias. Utilizar el objeto random para generar una posicion aleatoria de inicio
			// en la cadena sequence. Extraer la lectura de tamanho readLength e imprimirla en formato fastq.
			// Utilizar la cadena fixedQSStr para generar calidades fijas para el formato
			for (int i = 0; i < numReads; i++) {
			
				out.println("sequencia-" + (i + 1));
				char[] arregloSecuencias = sequence.substring(i, i + readLength).toCharArray();
				Set<Integer> a = new HashSet<>();
				//args[4]: Number of errors in a args[5] number of bp
				for (int j = 0; j < arregloSecuencias.length *(Double.valueOf(args[4])/Double.valueOf(args[5])); j++) 
				{
					int p = random.nextInt(arregloSecuencias.length);
				
					while (a.contains(p))
					{
						p = random.nextInt(arregloSecuencias.length);
						a.add(p);
					}
					char[] bases = new char[] { 'G', 'C', 'A', 'T' };
					char let = bases[random.nextInt(bases.length)];
					while (let == arregloSecuencias[p])
					{
						let = bases[random.nextInt(bases.length)];
					}
				}

				out.println(new String(arregloSecuencias));
				out.println(fixedQSStr);
				
			}
			out.close();
		}
	}
}
