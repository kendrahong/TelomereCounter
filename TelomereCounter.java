import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import org.apache.commons.lang3.StringUtils;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

public class TelomereCounter {
	/**
	 * 
	 * @author kendrahong
	 * Inputs: 
	 * args[0]: motif (e.g. TTAGGG)
	 * args[1]: telomeric repeat
	 * args[2]: FASTQ GZip textfile of telomeric reads 
	 * args[3]: directory for outputs
	 *
	 */
	public static void main(String[] args) throws Exception {
		long startTime;
		String motif = args[0];
		String oppositeMotif = findOpposite(motif);
		
		// Checks arguments
		if (!(new File(args[3]).isDirectory())) {
			throw new IllegalArgumentException(args[3] + "is not a directory. Raw output directory is required.");
		}
		
		// Reads data from each GZIP, counts telomeres, and prints results in the output file
		try (GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(args[2]))) {
			startTime = System.currentTimeMillis();
			BufferedReader br = new BufferedReader(new InputStreamReader(gzip));
			String line;
			long numberOfReads = 0;
			long elapsedTime = 0;
			Map<TelomereMarker, Long> telomereMap = new HashMap<TelomereMarker, Long>();
			for (int i = 0; i <= 5; i++) {
				for (int j = 0; j <= 5; j++) {
					telomereMap.put(new TelomereMarker(i, j, true), (long) 0);
					telomereMap.put(new TelomereMarker(i, j, false), (long) 0);
				}
			}

			long cytosineCount = 0;
			long guanineCount = 0;
			while ((line = br.readLine()) != null) {
				numberOfReads++;
				int firstSpace = line.indexOf("\t");
				int lastSpace = line.lastIndexOf("\t");
				String read = line.substring(firstSpace, lastSpace);
				int mid = read.length() / 2;

				// Counts telomeric sequences in each pair
				String[] pairs = { read.substring(0, mid), read.substring(mid, read.length()) };
				int pair1TelomereCount = 0;
				int pair2TelomereCount = 0;
				int readTTAGGGPair1 = StringUtils.countMatches(pairs[0], motif);
				int readTTAGGGPair2 = StringUtils.countMatches(pairs[1], motif);
				int readCCCTAAPair1 = StringUtils.countMatches(pairs[0], oppositeMotif);
				int readCCCTAAPair2 = StringUtils.countMatches(pairs[1], oppositeMotif);
				int readTTAGGG = readTTAGGGPair1 + readTTAGGGPair2;
				int readCCCTAA = readCCCTAAPair1 + readCCCTAAPair2;
				guanineCount += StringUtils.countMatches(read, "G");
				cytosineCount += StringUtils.countMatches(read, "C");

				// Count TTAGGG OR CCCTAA
				boolean cytosineORguanine = true;
				if (readTTAGGG == 0 && readCCCTAA == 0) {
					boolean guanine = read.contains("G");
					boolean cytosine = read.contains("C"); 
					if ((cytosine && !guanine) || guanine && !cytosine) {
						cytosineORguanine = false;
					}
				} else if (readTTAGGG == readCCCTAA) {
					pair1TelomereCount = readTTAGGGPair1;
					pair2TelomereCount = readTTAGGGPair2;
				} else if (readCCCTAA == 0 && readTTAGGG > 0) {
					pair1TelomereCount = readTTAGGGPair1;
					pair2TelomereCount = readTTAGGGPair2;
					cytosineORguanine = read.contains("C");
				} else if (readTTAGGG == 0 && readCCCTAA > 0){
					pair1TelomereCount = readCCCTAAPair1;
					pair2TelomereCount = readCCCTAAPair2;
					cytosineORguanine = read.contains("G");
				}
				for (TelomereMarker marker : telomereMap.keySet()) {
					if (marker.contains(pair1TelomereCount, pair2TelomereCount, cytosineORguanine)) {
						long val = telomereMap.get(marker);
						telomereMap.put(marker, val += 1);
						break;
					}
				}
			}

			// Timer
			long stopTime = System.currentTimeMillis();
			elapsedTime = stopTime - startTime;
	        long hours = TimeUnit.MILLISECONDS.toHours(elapsedTime);
	        elapsedTime -= TimeUnit.HOURS.toMillis(hours);
	        long minutes = TimeUnit.MILLISECONDS.toMinutes(elapsedTime);
	        elapsedTime -= TimeUnit.MINUTES.toMillis(minutes);
	        long seconds = TimeUnit.MILLISECONDS.toSeconds(elapsedTime);

			// Create new file in output folder and write results
			int index = args[2].lastIndexOf("/");
			String id = args[2].substring(index + 1).replaceFirst(".gz", "");

			// Sort the telomereMap by cytosine or guanine appearance
			long totalTelomere = 0;
			int telomericRepeat = Integer.parseInt(args[1]);
			int pairTelomericRepeat = telomericRepeat / 2;
			Map<TelomereMarker, Long> falseCytosineMap = new HashMap<TelomereMarker, Long>();
			Map<TelomereMarker, Long> trueCytosineMap = new HashMap<TelomereMarker, Long>();
			for (Map.Entry<TelomereMarker, Long> entry : telomereMap.entrySet()) {
				TelomereMarker markerEntry = entry.getKey();
				if (!markerEntry.getCytosine()) {
					falseCytosineMap.put(markerEntry, telomereMap.get(markerEntry));
					int pair1Count = markerEntry.getPair1Count();
					int pair2Count = markerEntry.getPair2Count();
					int totalCount = (pair1Count+pair2Count);
					if (pair1Count >= pairTelomericRepeat && pair2Count >= pairTelomericRepeat && totalCount >= telomericRepeat) {
						totalTelomere += (pair1Count+pair2Count) * entry.getValue();
					}
				} else {
					trueCytosineMap.put(markerEntry, telomereMap.get(markerEntry));
				}
			}

			// Sort the telomereMap by telomeric frequencies
			falseCytosineMap = sortByValues(falseCytosineMap);
			trueCytosineMap = sortByValues(trueCytosineMap);

			// Output the raw data
			File file = new File(args[3]  + id);
			try (Writer out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "utf-8"))) {
				out.write(String.format("%20s%n %20s%n %20s%n %20s%n %20s%n %20s%n %20s%n %20s%n", "ID: " + id.replace(".txt", ""),
					"NumberOfReads: " + numberOfReads, "ElapsedTime: " + hours + "hr" + minutes + "min" + seconds + "sec",
					"CytosineCount: " + cytosineCount, "GuanineCount: " + guanineCount, "Cytosine Probability: " + 
					(float)((float)(cytosineCount / 70)/(numberOfReads)), "Guanine Probability: " + 
					(float)((float)(guanineCount/70)/(numberOfReads)), "Telomeres with repeats >" + args[1] + ": " + totalTelomere));
				out.write(String.format("%20s %20s %20s %20s %20s%n", "Pair1Count", "Pair2Count", "Pair1Count+Pair2Count",
					"CytosineORGuanine", "TotalCount"));
				for (Map.Entry<TelomereMarker, Long> entry : falseCytosineMap.entrySet()) {
					TelomereMarker markerEntry = entry.getKey();
					out.write(String.format("%20d %20d %20d %20s %20d%n", markerEntry.getPair1Count(),
							markerEntry.getPair2Count(), markerEntry.getPair1Count() + markerEntry.getPair2Count(),
							markerEntry.getCytosine(), entry.getValue()));
				}
				for (Map.Entry<TelomereMarker, Long> entry : trueCytosineMap.entrySet()) {
					TelomereMarker markerEntry = entry.getKey();
					out.write(String.format("%20d %20d %20d %20s %20d%n", markerEntry.getPair1Count(),
							markerEntry.getPair2Count(), markerEntry.getPair1Count() + markerEntry.getPair2Count(),
							markerEntry.getCytosine(), entry.getValue()));
				}
				out.flush();
				out.close();
			} catch (IOException raw) {
				System.out.println("Unable to create raw output file inside " + args[3] + "\n" + raw.getMessage());
				System.exit(1);
			}
			System.out.println("Analysis of " + id + " is complete.");
		} catch (IOException e) {
			System.out.println("Error reading telomere GZIP files message for " + args[2] + "\n" + e.getMessage());
		}
	}

	// Sort the telomereMap by telomeric frequencies
	@SuppressWarnings("hiding")
	public static <TelomereMarker, Long extends Comparable<Long>> Map<TelomereMarker, Long> sortByValues(
			final Map<TelomereMarker, Long> map) {
		Comparator<TelomereMarker> valueComparator = new Comparator<TelomereMarker>() {
			public int compare(TelomereMarker m1, TelomereMarker m2) {
				int compare = map.get(m2).compareTo(map.get(m1));
				if (compare == 0)
					return 1;
				else
					return compare;
			}
		};
		Map<TelomereMarker, Long> sortedByValues = new TreeMap<TelomereMarker, Long>(valueComparator);
		sortedByValues.putAll(map);
		return sortedByValues;
	}
	
	public static String findOpposite(String motif) {
		String output = "";
		for (int i = motif.length() - 1; i >= 0; i--) {
			if (motif.charAt(i) == 'A') {
				output += "T";
			} else if (motif.charAt(i) == 'T') {
				output += "A";
			} else if (motif.charAt(i) == 'C') {
				output += "G";
			} else if (motif.charAt(i) == 'G') {
				output += "C";
			} else {
				throw new IllegalArgumentException("Invalid motif. Motifs can only contian ATCG"); 
			}
		}
		return output;	
	}
}
