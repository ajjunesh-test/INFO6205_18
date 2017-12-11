
package edu.neu.coe;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.IntStream;

public class KnapsackGAProblem {

	private boolean isMutated = false;
	private int crossoverCount = 0;
	private int cloneCount = 0;

	public int totalItems = 10;
	public int defaultPopulationSize = 10;
	public static final int sackCapacity = 50;
	public static final Double sackVolume = 3000.0;

	public double crossoverProbablity = 0.5;
	public double mutationProbabilty = 0.03;
	public int maxGen = 100;
	private int genCounter = 1;
	public double totalFitnessOfGeneration = 0;

	public ArrayList<Integer> chromosomeValueList;
	public ArrayList<Integer> chromosomeWeightList;
	public ArrayList<Double> chromosomeVolumeList;
	public ArrayList<Double> popFitness;
	private ArrayList<Double> bestFitnessList;
	private ArrayList<Double> avgFitnessList;
	public ArrayList<String> pop;
	public ArrayList<String> breededPop;
	private ArrayList<String> bestSolList;
	public Map<Integer, Double> popFitnessMap;

	public enum populationType {
		Normal, Breed
	}

	public static void main(String[] args) {

		KnapsackGAProblem kp = new KnapsackGAProblem();

	}

	public KnapsackGAProblem() {

		chromosomeValueList = new ArrayList<Integer>();
		chromosomeWeightList = new ArrayList<Integer>();
		chromosomeVolumeList = new ArrayList<Double>();
		popFitness = new ArrayList<Double>();
		bestFitnessList = new ArrayList<Double>();
		avgFitnessList = new ArrayList<Double>();
		pop = new ArrayList<String>();
		breededPop = new ArrayList<String>();
		bestSolList = new ArrayList<String>();
		popFitnessMap = new HashMap<Integer, Double>();

		this.generateInput();
		this.implementSolution();
		this.printFinalSolution();

	}

	public KnapsackGAProblem(String createEmptyObject) {
		
		chromosomeValueList = new ArrayList<Integer>();
		chromosomeWeightList = new ArrayList<Integer>();
		chromosomeVolumeList = new ArrayList<Double>();
		popFitness = new ArrayList<Double>();
		bestFitnessList = new ArrayList<Double>();
		avgFitnessList = new ArrayList<Double>();
		pop = new ArrayList<String>();
		breededPop = new ArrayList<String>();
		bestSolList = new ArrayList<String>();
		popFitnessMap = new HashMap<Integer, Double>();


	}

	public void implementSolution() {

		this.generateGeneList();

		System.out.println("\nGeneration 1:\n-------------------\nPopulation:");
		for (int i = 0; i < this.defaultPopulationSize; i++) {
			System.out.println((i) + " -- " + this.pop.get(i));
		}

		this.evaluatePopulation(populationType.Normal);
		Map<Integer,Double> newpopFitnessMap = sortByValue(popFitnessMap);

		System.out.println("\nFitness Score:\n---------------\n");
		for (Integer m : newpopFitnessMap.keySet()) {
			System.out.println((m) + " - " + newpopFitnessMap.get(m));
		}

		this.bestSolList.add(this.pop.get(this.bestSolPos()));

		System.out.println("\nBest solution:" + this.bestSolList.get(0));

		this.avgFitnessList.add(this.avgPopFitness());

		System.out.println("Average fitness: " + this.avgFitnessList.get(0));

		this.bestFitnessList.add(this.evaluateChromosome(this.pop.get(this.bestSolPos())));

		System.out.println("Best Fitness Score: " + this.bestFitnessList.get(0));

		if (this.maxGen > 1) {
			additionalGenerations();
		}

	}

	private void additionalGenerations() {

		for (int i = 1; i < this.maxGen; i++) {

			if ((this.maxGen > 5) && (i > 5)) {

				double a = this.avgFitnessList.get(i - 1);
				double b = this.avgFitnessList.get(i - 2);
				double c = this.avgFitnessList.get(i - 3);
				double d = this.avgFitnessList.get(i - 4);

				if (a == b && b == c && c == d) {
					maxGen = i;
					break;
				}
			}

			this.crossoverCount = 0;
			this.cloneCount = 0;
			this.isMutated = false;

			// Breed population normally
			/*
			 * for (int j = 0; j < this.populationSize / 2; j++) { this.breedPopulation(bo);
			 * }
			 */

			// Breed population parallelly
			IntStream.range(0, this.defaultPopulationSize / 2).parallel().forEach(counter -> {
				this.breedPopulation();
			});

			if (pop.size() != breededPop.size())
				this.breedPopulation();

			this.popFitness.clear();
			this.popFitnessMap.clear();

			this.evaluatePopulation(populationType.Breed);

			for (int k = 0; k < this.defaultPopulationSize; k++) {
				this.pop.set(k, this.breededPop.get(k));
			}

			System.out.println("\nGeneration " + (i + 1) + ":");
			System.out.println("***************\nPopulation:\n");
			for (int l = 0; l < this.defaultPopulationSize; l++) {
				System.out.println((l) + " - " + this.pop.get(l));
			}

			Map<Integer,Double> newpopFitnessMap = sortByValue(popFitnessMap);

			System.out.println("\nFitness Score:\n---------------\n");
			for (Integer m : newpopFitnessMap.keySet()) {
				System.out.println((m) + " - " + newpopFitnessMap.get(m));
			}

			this.breededPop.clear();

			this.bestSolList.add(this.pop.get(this.bestSolPos()));
			System.out.println("\nBest solution " + (i + 1) + ": " + this.bestSolList.get(i));

			this.avgFitnessList.add(this.avgPopFitness());
			System.out.println("Average fitness: " + this.avgFitnessList.get(i));

			this.bestFitnessList.add(this.evaluateChromosome(this.pop.get(this.bestSolPos())));

			System.out.println("Best Fitness: " + (i + 1) + ": " + this.bestFitnessList.get(i));

			System.out.println("Number of times Crossover Occured: " + this.crossoverCount);
			System.out.println("Number of times Cloning Occured: " + this.cloneCount);

			if (!this.isMutated)
				System.out.println("Mutation did not occur");
			else
				System.out.println("Mutation occured");

		}
	}

	private void printFinalSolution() {

		System.out.println("\nFinal items in knapsack: ");

		double bestFitness = 0;
		int bestGen = 0;

		for (int i = 0; i < this.maxGen - 1; i++) {
			if (this.bestFitnessList.get(i) > bestFitness) {
				bestFitness = this.bestFitnessList.get(i);
				bestGen = i;
			}
		}

		String bestSoln = this.bestSolList.get(bestGen);
		for (int i = 0; i < this.totalItems; i++) {
			if (bestSoln.substring(i, i+1).equalsIgnoreCase("1")) {
				System.out.print((i+1)+" ");
			}
		}

	}

	private void breedPopulation() {

		int chromosome1;
		int chromosome2;

		chromosome1 = selectChromosome();
		chromosome2 = selectChromosome();
		if (breededPop.size() == defaultPopulationSize)
			return;

		crossoverChromosomes(chromosome1, chromosome2);
		genCounter++;

	}

	public void mutateChromosome() {

		double rand = Math.random();
		if (rand <= mutationProbabilty) {

			this.isMutated = true;
			String mutChrom = "";
			String newMutChrom = "";
			int mutPt = 0;

			Random random = new Random();
			boolean selectFlag = random.nextBoolean();

			if (selectFlag) {
				mutChrom = breededPop.get(breededPop.size() - 1);
				mutPt = random.nextInt(totalItems);
				if (mutChrom.substring(mutPt, mutPt + 1).equals("1")) {
					newMutChrom = mutChrom.substring(0, mutPt) + "0" + mutChrom.substring(mutPt);

					breededPop.set(breededPop.size() - 1, newMutChrom);
				}
				if (mutChrom.substring(mutPt, mutPt + 1).equals("0")) {
					newMutChrom = mutChrom.substring(0, mutPt) + "1" + mutChrom.substring(mutPt);
					breededPop.set(breededPop.size() - 1, newMutChrom);
				}
				
			} else {
				mutChrom = breededPop.get(breededPop.size() - 2);
				mutPt = random.nextInt(totalItems);
				if (mutChrom.substring(mutPt, mutPt + 1).equals("1")) {
					newMutChrom = mutChrom.substring(0, mutPt) + "0" + mutChrom.substring(mutPt);
					breededPop.set(breededPop.size() - 1, newMutChrom);
				}
				if (mutChrom.substring(mutPt, mutPt + 1).equals("0")) {
					newMutChrom = mutChrom.substring(0, mutPt) + "1" + mutChrom.substring(mutPt);
					breededPop.set(breededPop.size() - 2, newMutChrom);
				}
				
			}
		}
	}

	public int selectChromosome() {

		double rand = Math.random() * totalFitnessOfGeneration;

		for (int i = 0; i < defaultPopulationSize; i++) {
			if (rand <= popFitnessMap.get(i)) {
				return i;
			}
			rand = rand - popFitnessMap.get(i);
		}

		return 0;
	}

	public void crossoverChromosomes(int chromosome1, int chromosome2) {

		String oldChr1 = pop.get(chromosome1);
		String oldChr2 = pop.get(chromosome2);
		StringBuffer newChr1 = new StringBuffer();
		StringBuffer newChr2 = new StringBuffer();

		double randCrossover = Math.random();
		if (randCrossover <= crossoverProbablity) {
			this.crossoverCount++;
			int count = 0;
			while (count < totalItems) {
				if (totalItems - count >= 4) {
					newChr1.append(oldChr1.charAt(count));
					newChr1.append(oldChr1.charAt(count + 1));
					newChr1.append(oldChr2.charAt(count + 2));
					newChr1.append(oldChr2.charAt(count + 3));

					newChr2.append(oldChr2.charAt(count));
					newChr2.append(oldChr2.charAt(count + 1));
					newChr2.append(oldChr1.charAt(count + 2));
					newChr2.append(oldChr1.charAt(count + 3));
					count += 4;
					continue;
				}
				if (totalItems - count >= 2) {
					newChr1.append(oldChr1.charAt(count));
					newChr1.append(oldChr2.charAt(count + 1));

					newChr2.append(oldChr2.charAt(count));
					newChr2.append(oldChr1.charAt(count + 1));
					count += 2;
					continue;
				}
				if (totalItems - count >= 1) {
					newChr1.append(oldChr1.charAt(count));
					newChr2.append(oldChr2.charAt(count));
					count += 1;
					continue;
				}
			}
			breededPop.add(String.valueOf(newChr1));
			breededPop.add(String.valueOf(newChr2));
		} else {
			this.cloneCount++;
			breededPop.add(pop.get(chromosome1));
			breededPop.add(pop.get(chromosome2));
		}

		mutateChromosome();
	}

	public int bestSolPos() {
		int bestPos = 0;
		double currentfitness = 0;
		double bestFitness = 0;
		for (int i = 0; i < defaultPopulationSize; i++) {
			currentfitness = evaluateChromosome(pop.get(i));
			if (currentfitness > bestFitness) {
				bestFitness = currentfitness;
				bestPos = i;
			}
		}
		return bestPos;
	}

	public double avgPopFitness() {
		double totFitness = 0;
		double avg = 0;
		for (int i = 0; i < defaultPopulationSize; i++) {
			totFitness = totFitness + popFitnessMap.get(i);
		}
		avg = totFitness / defaultPopulationSize;
		return avg;
	}

	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> popFitnessMap) {
		List<Map.Entry<K, V>> list = new ArrayList<Map.Entry<K, V>>(popFitnessMap.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<K, V>>() {
			public int compare(Map.Entry<K, V> o1, Map.Entry<K, V> o2) {
				return (o1.getValue()).compareTo(o2.getValue());
			}
		});

		Map<K, V> result = new LinkedHashMap<K, V>();
		for (Map.Entry<K, V> entry : list) {
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}

	public void evaluatePopulation(populationType type) {
		totalFitnessOfGeneration = 0;
		for (int i = 0; i < defaultPopulationSize; i++) {
			double tFitness;
			if (type.equals(populationType.Normal))
				tFitness = evaluateChromosome(pop.get(i));
			else
				tFitness = evaluateChromosome(breededPop.get(i));
			popFitness.add(tFitness);
			popFitnessMap.put(i, tFitness);
			totalFitnessOfGeneration = totalFitnessOfGeneration + tFitness;
		}
	}

	public double evaluateChromosome(String chromosome) {
		double chromosomeWt = 0;
		double chromosomeVal = 0;
		double chromosomeVol = 0;
		double fitness = 0;
		char c = '0';

		for (int j = 0; j < totalItems; j++) {
			c = chromosome.charAt(j);
			if (c == '1') {
				chromosomeWt = chromosomeWt + chromosomeWeightList.get(j);
				chromosomeVal = chromosomeVal + chromosomeValueList.get(j);
				chromosomeVol = chromosomeVol + chromosomeVolumeList.get(j);

			}
		}
		if (sackCapacity >= chromosomeWt && sackVolume >= chromosomeVol)
			fitness = chromosomeVal;
		return fitness;
	}

	private void generateGeneList() {
		for (int i = 0; i < defaultPopulationSize; i++) {
			pop.add(makeChromosome());
		}
	}

	public String makeChromosome() {

		StringBuilder chhromosome = new StringBuilder(totalItems);
		char c;
		for (int i = 0; i < totalItems; i++) {
			c = '0';
			double rnd = Math.random();
			if (rnd > 0.5)
				c = '1';
			chhromosome.append(c);
		}
		return chhromosome.toString();
	}

	private void generateInput() {

		System.out.println("The number of items: " + totalItems);

		Random rand = new Random();
		for (int i = 0; i < totalItems; i++) {
			chromosomeValueList.add(i + 1);
			Double volume;
			int weight = rand.nextInt(sackCapacity / 2);
			System.out.println("the weight of item " + (i + 1) + ": " + weight);
			if (weight != 0)
				chromosomeWeightList.add(weight);
			else
				chromosomeWeightList.add(1);

			if (weight != 0)
				volume = weight * rand.nextDouble() * 100;
			else
				volume = 1 * rand.nextDouble() * 100;
			System.out.println("the volume of item " + (i + 1) + ": " + volume);
			chromosomeVolumeList.add(volume);

		}

		// Capacity of knapsack
		System.out.println("The knapsack capacity: " + sackCapacity);

		// Volume of knapsack
		System.out.println("The knapsack volume: " + sackVolume);

		// Population size
		System.out.println("The population size: " + defaultPopulationSize);

		// Maximum number of generations
		System.out.println("The maximum number of generations: " + maxGen);

		System.out.println("The crossover probability: " + crossoverProbablity);
		System.out.println("The mutation probability: " + mutationProbabilty);

	}
}