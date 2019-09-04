package uniandes.algorithms.readsanalyzer;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import ngsep.sequences.RawRead;
/**
 * Stores abundances information on a list of subsequences of a fixed length k (k-mers)
 * @author Jorge Duitama
 */
public class KmersTable implements RawReadProcessor {

	private int kmerLenght=0;
	private Map<String, Integer> kmerCount;

	/**
	 * Creates a new table with the given k-mer size
	 * @param kmerSize length of k-mers stored in this table
	 */
	public KmersTable(int kmerSize) {
		// TODO: Implementar metodo

		this.kmerLenght=kmerSize;
		kmerCount= new HashMap<String, Integer>();
	}

	/**
	 * Identifies k-mers in the given read
	 * @param read object to extract new k-mers
	 */
	public void processRead(RawRead read) {
		String sequence = read.getSequenceString();
		// TODO Implementar metodo. Calcular todos los k-mers del tamanho dado en la constructora y actualizar la abundancia de cada k-mer

		for (int i = 0; i <= sequence.length() - kmerLenght; i++) {
			String kmer = sequence.substring(i, i + kmerLenght);
			if (kmerCount.containsKey(kmer))
			{	kmerCount.compute(kmer, (key, j) -> j + 1);
			}
			else
				kmerCount.put(kmer, 1);
		}
	}

	/**
	 * List with the different k-mers found up to this point
	 * @return Set<String> set of k-mers
	 */
	public Set<String> getDistinctKmers() {
		// TODO Implementar metodo
		return kmerCount.keySet();
	}

	/**
	 * Calculates the current abundance of the given k-mer 
	 * @param kmer sequence of length k
	 * @return int times that the given k-mer have been extracted from given reads
	 */
	public int getAbundance(String kmer) {
		// TODO Implementar metodo
		return kmerCount.get(kmer);
	}

	/**
	 * Calculates the distribution of abundances
	 * @return int [] array where the indexes are abundances and the values are the number of k-mers
	 * observed as many times as the corresponding array index. Position zero should be equal to zero
	 */
	public int[] calculateAbundancesDistribution() {
		// TODO Implementar metodo

		int[] distributionAbundances = new int[kmerCount.size()];
		kmerCount.values().stream().forEach((c) -> distributionAbundances[c]++);
		return distributionAbundances;
	}
}
