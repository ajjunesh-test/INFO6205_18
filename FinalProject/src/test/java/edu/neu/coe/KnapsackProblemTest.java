package edu.neu.coe;
import org.junit.Test;

import edu.neu.coe.KnapsackGAProblem;
import junit.framework.TestCase;
public class KnapsackProblemTest extends TestCase {
	
	public KnapsackGAProblem kp;
	
	@Override
	protected void setUp() {
		kp = new KnapsackGAProblem("test");
		addPopulation();
		setItemValueAndWeight();
		
		
	}
	
	
	public void addPopulation() {
		kp.pop.add("1111011000");
		kp.pop.add("1111011011");
		kp.pop.add("0000010100");
		kp.pop.add("0101100111");
		kp.pop.add("0110100101");
		kp.pop.add("0100110000");
		kp.pop.add("1001110001");
		kp.pop.add("1011011110");
		kp.pop.add("1101101100");
		kp.pop.add("0101101100");
		
	}
	
	public void setItemValueAndWeight() {
		for(int i=0;i<kp.totalItems;i++) {
			kp.chromosomeValueList.add(i+1);
			kp.chromosomeWeightList.add(i+1);
			kp.chromosomeVolumeList.add((i+1)*1.0*25);
		}
	}
	
	
	@Test
	public void testMakeGene() {
		String newGene = kp.makeChromosome();
		assertEquals(10, newGene.length());
		System.out.println("testMakeGene() run successfull");
	}
	
	
	@Test
	public void testEvaluateGene() {
		
		String gene1 = "1111011000";
		String gene2 = "1111011011";
		
		Double fitness = kp.evaluateChromosome(gene1);
		Double fitness2 = kp.evaluateChromosome(gene2);
		
		assertEquals(23.0, fitness);
		assertEquals(42.0, fitness2);
		System.out.println("testEvaluateGene() run successfull");
		
	}
	
	@Test
	public void testEvaluatePopulation() {
		kp.evaluatePopulation(KnapsackGAProblem.populationType.Normal);
		assertEquals(275.0, kp.totalFitnessOfGeneration);
		System.out.println("testEvaluatePopulation() run successfull");
	}
	
	
	@Test
	public void testAvgFitness() {
		kp.evaluatePopulation(KnapsackGAProblem.populationType.Normal);
		assertEquals(27.5, kp.avgPopFitness());
		System.out.println("testAvgFitness() run successfull");
	}
	
	@Test
	public void testGetBestSolutionPosition() {
		assertEquals(1, kp.bestSolPos());
		System.out.println("testGetBestSolutionPosition() run successfull");
		
	}
	
	@Test
	public void testSelectChromosome() {
		boolean rangeFlag = false;
		kp.evaluatePopulation(KnapsackGAProblem.populationType.Normal);
		int selectedChromosome = kp.selectChromosome();
		
		if(0 < selectedChromosome && selectedChromosome < kp.defaultPopulationSize)
			rangeFlag = true;
		
		assertEquals(true, rangeFlag);
		System.out.println("testSelectChromosome() runs successfull");
		
	}
	
	@Test
	public void testCrossover() {
		int gene1 = 0;
		int gene2 = 1;
		
		kp.crossoverProbablity = 1.0;
		kp.mutationProbabilty = 0.0;
		kp.crossoverChromosomes(gene1, gene2);
		
		assertEquals("1111011001", kp.breededPop.get(0));
		assertEquals("1111011010", kp.breededPop.get(1));
			
		System.out.println("testCrossover() runs successfull");
		
	}
	
	@Test
	public void testmutateChromosome() {
		
		kp.mutationProbabilty = 1.0;
		kp.breededPop.clear();
		kp.totalItems = 2;
		kp.defaultPopulationSize = 2;
		kp.breededPop.add("00");
		kp.breededPop.add("11");
		
		kp.mutateChromosome();

		assertTrue(kp.breededPop.get(0).equals("00") || 
				kp.breededPop.get(0).equals("010") ||
				kp.breededPop.get(0).equals("100"));
		
		assertTrue(kp.breededPop.get(1).equals("011") || 
				kp.breededPop.get(1).equals("11") ||
				kp.breededPop.get(1).equals("11"));
		
		//cases for 00, 11
		//case 1 -> 00,011
		//case 2 -> 010,11
		//case 3 -> 00,101
		//case 4 -> 100,11

		System.out.println("testmutateChromosome() runs successfull");
		
		
	}

}
