package uniandes.algorithms.readsanalyzer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.Map.Entry;

import htsjdk.samtools.util.RuntimeEOFException;
import ngsep.math.Distribution;
import ngsep.sequences.RawRead;

/**
 * Represents an overlap graph for a set of reads taken from a sequence to assemble
 * @author Jorge Duitama, Christian Chavarro, Daniel Bautista
 *
 */
public class OverlapGraph implements RawReadProcessor {

	private int minOverlap;
	private Map<String,Integer> readCounts = new HashMap<>();
	private Map<String,ArrayList<ReadOverlap>> overlaps = new HashMap<>();

	/**
	 * Creates a new overlap graph with the given minimum overlap
	 * @param minOverlap Minimum overlap
	 */
	public OverlapGraph(int minOverlap) {
		this.minOverlap = minOverlap;
	}

	/**
	 * Adds a new read to the overlap graph
	 * @param read object with the new read
	 */
	public void processRead(RawRead read) {
		String sequence = read.getSequenceString();
		//TODO: Paso 1. Agregar la secuencia al mapa de conteos si no existe.
		//Si ya existe, solo se le suma 1 a su conteo correspondiente y no se deben ejecutar los pasos 2 y 3 
		if (!readCounts.containsKey(sequence))
		{
			//Se agrega en la primera posici�n del mapa
			readCounts.put(sequence, 1);
			//La lista de sobrelapes en las secuencias es creada. Puede ser un Array, se deja el tipo generico "list"
			ArrayList<ReadOverlap> overlapSequence = new ArrayList<ReadOverlap>();
			//Se crea la variable SequenceLenght para no solicitar dos veces el numero de la longitud de la secuencia. Se consume m�s espacio en memoria
			// Pero se reduce el tiempo de computo
			ReadOverlap overlapLecture = null;
			int sequenceLenght = sequence.length();
			String suf =sequence.substring(sequenceLenght-minOverlap, sequenceLenght-1);
			String pre = "";

			for (String sequencePrefix:readCounts.keySet())
			{
				pre = sequencePrefix.substring(0, suf.length()-1);
				if (pre.equals(suf)){
					overlapLecture= new ReadOverlap(sequence, sequencePrefix, suf.length());
					overlapSequence.add(overlapLecture);
				}
			}
			overlaps.put(sequence, overlapSequence);

		}
		else
		{
			readCounts.put(sequence, readCounts.get(sequence)+1);
		}

		//TODO: Paso 2. Actualizar el mapa de sobrelapes con los sobrelapes en los que la secuencia nueva sea predecesora de una secuencia existente
		//2.1 Crear un ArrayList para guardar las secuencias que tengan como prefijo un sufijo de la nueva secuencia
		//2.2 Recorrer las secuencias existentes para llenar este ArrayList creando los nuevos sobrelapes que se encuentren.
		//2.3 Después del recorrido para llenar la lista, agregar la nueva secuencia con su lista de sucesores al mapa de sobrelapes 

		//TODO: Paso 3. Actualizar el mapa de sobrelapes con los sobrelapes en los que la secuencia nueva sea sucesora de una secuencia existente
		// Recorrer el mapa de sobrelapes. Para cada secuencia existente que tenga como sufijo un prefijo de la nueva secuencia
		//se agrega un nuevo sobrelape a la lista de sobrelapes de la secuencia existente


	}
	/**
	 * Returns the length of the maximum overlap between a suffix of sequence 1 and a prefix of sequence 2
	 * @param sequence1 Sequence to evaluate suffixes
	 * @param sequence2 Sequence to evaluate prefixes
	 * @return int Maximum overlap between a prefix of sequence2 and a suffix of sequence 1
	 */
	private int getOverlapLength(String sequence1, String sequence2) {
		// TODO Implementar metodo
		int lenght =0;
		String prefix="";
		boolean find =false;
		int s2Lenght = sequence2.length();
		for(int i=0; i<s2Lenght&&!find;i++)
		{
			prefix=sequence2.substring(0,s2Lenght-1-i);
			if(sequence1.endsWith(prefix))
			{
				lenght=prefix.length();
				find=true;
			}
		}
		return lenght;
	}



	/**
	 * Returns a set of the sequences that have been added to this graph 
	 * @return Set<String> of the different sequences
	 */
	public Set<String> getDistinctSequences() {
		//TODO: Implementar metodo
		return overlaps.keySet();
	}

	/**
	 * Calculates the abundance of the given sequence
	 * @param sequence to search
	 * @return int Times that the given sequence has been added to this graph
	 */
	public int getSequenceAbundance(String sequence) {
		//TODO: Implementar metodo
		return readCounts.get(sequence);
	}

	/**
	 * Calculates the distribution of abundances
	 * @return int [] array where the indexes are abundances and the values are the number of sequences
	 * observed as many times as the corresponding array index. Position zero should be equal to zero
	 */
	public int[] calculateAbundancesDistribution() {
		//TODO: Implementar metodo
		int[] distribution = new int[readCounts.size()];
		int frequency=0;
		for(String sequence:readCounts.keySet())
		{
			frequency=getSequenceAbundance(sequence);
			distribution[frequency]=distribution[frequency]+1;
		}
		return distribution;
	}
	/**
	 * Calculates the distribution of number of successors
	 * @return int [] array where the indexes are number of successors and the values are the number of 
	 * sequences having as many successors as the corresponding array index.
	 */
	public int[] calculateOverlapDistribution() {
		// TODO: Implementar metodo
		int[] distribution = new int[overlaps.size()];
		int frequency=0;
		for(String sequence:overlaps.keySet())
		{
			frequency=overlaps.get(sequence).size();
			distribution[frequency]=distribution[frequency]+1;
		}
		return distribution;
	}
	/**
	 * Predicts the leftmost sequence of the final assembly for this overlap graph
	 * @return String Source sequence for the layout path that will be the left most subsequence in the assembly
	 */
	public String getSourceSequence () {
		// TODO Implementar metodo recorriendo las secuencias existentes y buscando una secuencia que no tenga predecesores
		Map<String, Integer> sequencePredecesors = new HashMap<String, Integer>();
		for (ArrayList<ReadOverlap> listReadOverlaps : overlaps.values()) 
		{
			for (ReadOverlap overlaps : listReadOverlaps) 
			{
				String act = overlaps.getDestSequence();
				if (sequencePredecesors.containsKey(act)) 
				{
					sequencePredecesors.compute(act, (key, i) -> i + 1);
				} 
				else
				{	
					sequencePredecesors.put(act, 1);
				}
			}
		}
		String r=sequencePredecesors.entrySet().stream().min((a, b) -> a.getValue() - b.getValue()).get().getKey();
System.out.println("Entra al sequence source" + r);
		return r;
	}

	/**
	 * Calculates a layout path for this overlap graph
	 * @return ArrayList<ReadOverlap> List of adjacent overlaps. The destination sequence of the overlap in 
	 * position i must be the source sequence of the overlap in position i+1. 
	 */
	public ArrayList<ReadOverlap> getLayoutPath() {
		ArrayList<ReadOverlap> layout = new ArrayList<>();
		HashSet<String> visitedSequences = new HashSet<>(); 
		// TODO Implementar metodo. Comenzar por la secuencia fuente que calcula el método anterior
		// Luego, hacer un ciclo en el que en cada paso se busca la secuencia no visitada que tenga mayor sobrelape con la secuencia actual.
		// Agregar el sobrelape a la lista de respuesta y la secuencia destino al conjunto de secuencias visitadas. Parar cuando no se encuentre una secuencia nueva
		String sequence=getSourceSequence();
		while (true) 
		{
			visitedSequences.add(sequence);
			Optional<ReadOverlap> next = overlaps.get(sequence).stream().filter((a) -> !visitedSequences.contains(a.getDestSequence())).max((a, b) -> a.getOverlap() - b.getOverlap());

			if (next.isPresent()) 
			{
				sequence = next.get().getDestSequence();
				layout.add(next.get());
			}
			else
			{
				break;
			}
		}
		return layout;
	}
	/**
	 * Predicts an assembly consistent with this overlap graph
	 * @return String assembly explaining the reads and the overlaps in this graph
	 */
	public String getAssembly () {
		ArrayList<ReadOverlap> layout = getLayoutPath();
		StringBuilder assembly = new StringBuilder();
		// TODO Recorrer el layout y ensamblar la secuencia agregando al objeto assembly las bases adicionales que aporta la región de cada secuencia destino que está a la derecha del sobrelape 
		assembly.append(layout.get(0).getSourceSequence());

		for (ReadOverlap response:layout)
		{
			assembly.append(response.getDestSequence().substring(response.getOverlap()));
		}

		return assembly.toString();
	}


}
