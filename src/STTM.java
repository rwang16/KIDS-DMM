
import eval.TopicQualityEval;
import models.*;


import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;

import utility.CmdArgs;
import eval.ClusteringEval;
import eval.ClassificationEval;

/**
 * STTM: A Java package for the short text topic models including DMM, BTM, WNTM, PTM, SATM, LDA, LFDMM, LFLDA, etc.
 *
 *
 * @author: Jipeng Qiang
 *
 * @version: 1.0
 *
 */
public class STTM
{
    public static void main(String[] args)
    {

        CmdArgs cmdArgs = new CmdArgs();
        CmdLineParser parser = new CmdLineParser(cmdArgs);
        try {

            parser.parseArgument(args);

            if (cmdArgs.model.equals("KIDSDMM")) {
                KIDSDMM bdc = new KIDSDMM(cmdArgs.corpus, cmdArgs.output_dir, cmdArgs.schema, cmdArgs.threshold,
                        cmdArgs.filterSize, cmdArgs.ntopics, cmdArgs.alpha, cmdArgs.beta,
                        cmdArgs.lambda, cmdArgs.niters,cmdArgs.nBurnIn, cmdArgs.twords, cmdArgs.expModelName);
                bdc.inference();
            } else if (cmdArgs.model.equals("ClusteringEval")) {
                ClusteringEval.evaluate(cmdArgs.labelFile, cmdArgs.dir,
                        cmdArgs.prob);
            }else if (cmdArgs.model.equals("ClassificationEval")) {
                ClassificationEval.evaluate(cmdArgs.labelFile, cmdArgs.dir,
                        cmdArgs.prob);
            }else if (cmdArgs.model.equals("TopicQualityEval")) {
                TopicQualityEval topicQualityEval = new TopicQualityEval(cmdArgs.ntopics, cmdArgs.topWordsPath,
                        cmdArgs.vocabFile, cmdArgs.wikiDir, cmdArgs.topTC, cmdArgs.topicCoherEval);
                topicQualityEval.computeTopicQualityEval();
            }
            else {
                System.out.println("Error: Option \"-model\" must get \"KIDSDMM\" or " +
                        "\"ClusteringEval\" or \"ClassificationEval\" or \"TopicQualityEval\"");
                System.out
                        .println("\tKIDSDMM: Infer topics using KIDSDMM");
                System.out
                        .println("\tClusteringEval: Specify the document clustering evaluation");
                System.out
                        .println("\tClassificationEval: Specify the document classification evaluation");
                System.out
                        .println("\tTopicQualityEval: Specify the topic quality evaluation");
                help(parser);
                return;
            }
        }
        catch (CmdLineException cle) {
            System.out.println("Error: " + cle.getMessage());
            help(parser);
            return;
        }
        catch (Exception e) {
            System.out.println("Error: " + e.getMessage());
            e.printStackTrace();
            return;
        }

        System.out.println("end!!!!!!!");
    }

    public static void help(CmdLineParser parser)
    {
        System.out
                .println("java -jar jSTTM.jar [options ...] [arguments...]");
        parser.printUsage(System.out);
    }
}
