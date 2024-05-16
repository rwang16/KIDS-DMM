//package utility;
//
//import java.util.*;
//
//public class TermWeighting {
//    public double[] termWeight;
//    public List<List<Integer>> corpus;
//    public int vocabularySize;
//    public int numTopic;
//    public int numWordsInCorpus;
//
//    public double[][] phi;
//    public int[][] topicTopWords;
//
//    public int lambda;
//
//    public double[] logTermWeight() {
//        double[] logWeight = new double[vocabularySize];
//
//        int[] wordCount = new int[vocabularySize];
//        for (List<Integer> doc: corpus) {
//            for (int wordId: doc) {
//                wordCount[wordId] += 1;
//            }
//        }
//        for (int i = 0; i < vocabularySize; i++) {
//            if (wordCount[i] == 1) {
//                logWeight[i] = -Math.log(1.0 / vocabularySize);
//            } else {
//                logWeight[i] = -Math.log(1.0 * wordCount[i] / numWordsInCorpus);
//            }
//        }
//        return logWeight;
//    }
//
//    public double[] pmiTermWeight() {
//        double[] pmiWeight = new double[vocabularySize];
//
//        int[] wordCount = new int[vocabularySize];
//        Map<Integer, Set<Integer>> wordOccDocCount = new HashMap<>();
//        for (int docId = 0; docId < corpus.size(); docId++) {
//            for (int wordId: corpus.get(docId)) {
//                Set<Integer> docSet = new HashSet<>();
//                wordCount[wordId] += 1;
//                if (wordOccDocCount.containsKey(wordId)) {
//                    docSet = wordOccDocCount.get(wordId);
//                }
//                docSet.add(docId);
//                wordOccDocCount.put(wordId, docSet);
//            }
//        }
//        for (int i = 0; i < vocabularySize; i++) {
//            pmiWeight[i] = -Math.log(1.0 * wordOccDocCount.get(i).size()
//                    / corpus.size() / wordCount[i] * numWordsInCorpus);
//        }
//        return pmiWeight;
//    }
//
//    public double[] bdcTermWeight() {
//        double[] bdcWeight = new double[vocabularySize];
//
//        double[] phiSumTopic = new double[vocabularySize];
//        for (int i = 0; i < vocabularySize; i++) {
//            double sum = 0;
//            for (int tIndex = 0; tIndex < numTopic; tIndex++) {
//                sum += phi[tIndex][i];
//            }
//            phiSumTopic[i] = sum;
//        }
//        for (int i = 0; i < vocabularySize; i++) {
//            double sum = 0.0;
//            for (int j = 0; j < numTopic; j++) {
//                double pro = phi[j][i] / phiSumTopic[i];
//                sum += pro * Math.log(pro);
//                if (Double.isNaN(sum)) {
//                    System.out.println("***************");
//                }
//            }
//            bdcWeight[i] = 1 + sum / Math.log(numTopic);
//        }
//        return bdcWeight;
//    }
//
//    public double[] ewTermWeight() {
//        double[] ewWeight = new double[vocabularySize];
//
//        int[][] wordCoOccur = new int[vocabularySize][vocabularySize];
//        int[] wordPairCount = new int[vocabularySize];
//        for (List<Integer> doc: corpus){
//            for (Integer i: doc) {
//                for (Integer j: doc) {
//                    if (!i.equals(j)) {
//                        wordCoOccur[i][j] += 1;
//                        wordCoOccur[j][i] += 1;
//                        wordPairCount[i] += 1;
//                        wordPairCount[j] += 1;
//                    }
//                }
//            }
//        }
//        for (int i = 0; i < vocabularySize; i++) {
//            double sum = 0.0;
//            for (int j = 0; j < vocabularySize; j++) {
//                if (i != j && wordCoOccur[i][j] != 0.0) {
//                    double pro = 1.0 * wordCoOccur[i][j] / wordPairCount[i];
//                    sum += pro * Math.log(pro);
//                }
//            }
//            if (sum == 0.0) {
//                ewWeight[i] = 1.0;
//            } else{
//                ewWeight[i] = (-sum) / Math.log(wordPairCount[i]);
//            }
//        }
//
//        return ewWeight;
//    }
//
//    public double[] mTFIDFTermWeight() {
//        double[] mTFIDFWeight = new double[vocabularySize];
//        double[] wIDF = new double[vocabularySize];
//        for (int i = 0; i < wIDF.length; i++) {
//            wIDF[i] = Math.log(1. * corpus.size() / wDF[i]) + 1;
//        }
//
//        for (int docId = 0; docId < Corpus.size(); docId++) {
//            int[] doc = Corpus.get(docId);
//            Map<Integer, Integer> wordOccurs = new HashMap<>();
//            for (int wordId: doc) {
//                if (wordOccurs.containsKey(wordId)) {
//                    int times = wordOccurs.get(wordId);
//                    times += 1;
//                    wordOccurs.put(wordId, times);
//                } else {
//                    wordOccurs.put(wordId, 1);
//                }
//            }
//            int sum = 0;
//            for (int wordId: wordOccurs.keySet()) {
//                sum += Math.pow(wordOccurs.get(wordId), 2);
//            }
//            for (int i = 0; i < doc.length; i++) {
//                int wordId = doc[i];
//                mTFIDF[wordId] += wordOccurs.get(wordId) * Math.log(Math.sqrt(vocabularySize) / termCount[wordId])
//                        / Math.log(sum * Math.pow(wordOccurs.size(),2) / Math.sqrt(vocabularySize));
//            }
//        }
//
//        for (int i = 0; i < vocabularySize; i++) {
//            mTFIDF[i] = mTFIDF[i] * wIDF[i];
//        }
//        return mTFIDFWeight;
//    }
//
//    public double[] iwpTermWeight() {
//        double[] iwpWeight = new double[vocabularySize];
//        // Count the num of word occurrence
//        Map<Integer, Set<Integer>> cnt = new HashMap<>();
//        for (int i = 0; i < topicTopWords.length; i++) {
//            for (int j = 0; j < topicTopWords[0].length; j++) {
//                int wordId = topicTopWords[i][j];
//                Set<Integer> topicSet = new HashSet<>();
//                if (cnt.containsKey(wordId)){
//                    topicSet = cnt.get(wordId);
//                }
//                topicSet.add(i);
//                cnt.put(wordId, topicSet);
//            }
//        }
//        Set<Integer> wordIds = cnt.keySet();
//        // Compute the mean of TU of top word
//        for (int wordId: wordIds) {
//            Set<Integer> topicSet = cnt.get(wordId);
//            iwpWeight[wordId] = Math.exp(1.0 * (lambda - topicSet.size()) / numTopic);
//
//        }
//        return iwpWeight;
//    }
//}
