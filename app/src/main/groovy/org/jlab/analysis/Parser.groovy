package org.jlab.analysis;

// Groovy Imports
import groovy.transform.CompileStatic;

// Java Imports
import java.io.*;
import java.util.*;
import java.text.*;
import org.jlab.clas.physics.LorentzVector;

/**
* Encapsulates most of the command line parsing and messaging for the Analysis class.
*
* @version 1.0
* @author  Matthew McEneaney
*/

@CompileStatic
public class Parser {

    /**
    * Constructor stub
    */
    public Parser() {}

    /**
    * Print out help message.  Pretty straightforward right?
    */
    protected boolean help() {

        System.out.println("An exception occurred: Use the --help option to see proper usage.");
        System.out.println("------------------------------------------------------------------------------");
        return false;
    }

    /**
    * Print out author, contact, version, etc.
    *
    * @return Boolean flag
    */
    protected boolean version() {

        System.out.println(" Author:  Matthew McEneaney ");
        System.out.println(" Contact: matthew.mceneaney@duke.edu ");
        System.out.println(" Version: 1.0");
        System.out.println(" Groovy:  " + GroovySystem.getVersion());
        System.out.println("------------------------------------------------------------------------------");
        return false;
    }

    /**
    * Print out available kinematic variables for an Analysis
    *
    * @param analysis
    *
    * @return Boolean flag
    */
    protected boolean kinematics(Analysis analysis) {

        System.out.println(" Available Kinematics: ");
        String spacing = "  "; for (String n : analysis.getKinematics().keySet()) {spacing += n+" ";}
        System.out.println(spacing);
        System.out.println("------------------------------------------------------------------------------");
        return false;
    }

    /**
    * Print out available command line options with descriptions.
    *
    * @return Boolean flag
    */
    protected boolean options() {

        System.out.println(" Usage : Analysis [input file or directory] [options] ");
        System.out.println("\n * General Options :");
        System.out.println("\t-vtx           : Add particle vertices to tree");
        System.out.println("\t-ang           : Add particles' polar/azimuthal angles to tree");
        System.out.println("\t-q             : Identify particles by charge");
        System.out.println("\t-e             : Do not require scattered electron in FT");
        System.out.println("\t-s             : Use strict pid to mass assignment for kinematics");
        System.out.println("\t-rn            : Include run # in NTuple");
        System.out.println("\t-en            : Include event # in NTuple");
        System.out.println("\t-lk            : Include extra kinematics for Lambda analysis");
        System.out.println("\t-ak            : Include Affinity partonic kinematics for MC matched-events");
        System.out.println("\t-ik            : Include extra kinematics for individual particles");
        System.out.println("\t-be   [float]  : Beam energy (GeV) for kinematics   (def: 10.6)");
        System.out.println("\t-tm   [float]  : Target mass (GeV) for kinematics   (def: 0.9383)");
        System.out.println("\t-tpid [int]    : Target Lund pid for kinematics     (def: 2212)");
        System.out.println("\t-tspin [float..] : Target spin vector components for kinematics (def: 0 0 1)");
        System.out.println("\t-tspin_sign [int] : Target spin vector sign (def: 1)");
        System.out.println("\t-qa            : Use clasqaDB check");
        System.out.println("\t-set_run [int] : Set run # in NTuple");
        System.out.println("\t-qm   [string] : Specify clasqaDB check method (def: OkForAsymmetry)");
        System.out.println("\t-fc   [int]    : Use fiducial volume cuts (0-loose, 1-med, 2-tight)");
        System.out.println("\t-momc [string, int, int] : Use momentum corrections. Parameters:");
        System.out.println("\t                           * dataset (Fall2018 or Spring2019)");
        System.out.println("\t                           * pass version (1,2,3,...)");
        System.out.println("\t                           * outbending (0 : false, 1 : true)");
        System.out.println("\t-mcsmear [string, double, double, double, bool] : Use MC smearing. Parameters:");
        System.out.println("\t                           * path to JSON file containing:");
        System.out.println("\t                               - 'mombinlims' : Map of PIDs to bin ids to momentum bin limits");
        System.out.println("\t                               - 'mom'        : Map of PIDs to bin ids to delta p / p means and widths");
        System.out.println("\t                               - 'theta'      : Map of PIDs to bin ids to delta theta means and widths");
        System.out.println("\t                               - 'phi'        : Map of PIDs to bin ids to delta phi means and widths");
        System.out.println("\t                           * Momentum additional data smearing fraction");
        System.out.println("\t                           * Theta additional data smearing fraction");
        System.out.println("\t                           * Phi additional data smearing fraction");
        System.out.println("\t                           * Option to offset by means (mus) of smearing distributions");
        System.out.println("\t-xF   [float]  : Minimum xF cut for events          (def: none)");
        System.out.println("\t-m    [float]  : Maximum mass cut for events        (def: none)");
        System.out.println("\t-'[var]>(<)[float]' : Cut kinematic variable");
        System.out.println("\n * PID Tag options");
        System.out.println("\t-tag  [int...]       : Require pid tag(s) in event");
        System.out.println("\t-filter [int,int...] : Filter events by pid counts <=count if filter>0");
        System.out.println("\t                       or >count if filter<0 (delimiter ':')");
        System.out.println("\n * ML Client Options");
        System.out.println("\t-mlclient  [hostname, port] : Connect to ML model server (default localhost,5000)");
        System.out.println("\t-mlclient_banks  [input bank names] : ML model server input bank names (delimiter ',', default '')");
        System.out.println("\t-mlclient_nscores  [int] : ML model server number of output scores (default 2)");
        System.out.println("\n * Channel/Bank Options");
        System.out.println("\t-ch   [int...] : Specify Lund pids to look for in REC::Particle bank");
        System.out.println("\t                  (delimiter ':', default 2212:-211).");
        System.out.println("\t                  Denote groups for kinematics by just separating with commas,");
        System.out.println("\t                  e.g., 22:2212,-211:321 groups 2212 and -211.");
        System.out.println("\t-mch  [int...] : Specify Lund pids to look for in MC::Lund bank");
        System.out.println("\t-pch  [int...] : Specify parent pids to look for in MC::Lund bank");
        System.out.println("\t-ma            : Require matching decay in MC::Lund bank");
        //TODO: Add mixing options and sector cut options
        //TODO: Option to include parents in tree or not -> Methods for reading all momenta and vertices or not?  And -> Finish mcmatch methods match just copies decay to mcchan
        // System.out.println("\t-mc           : Look for decay in MC::Lund bank");
        System.out.println("\n * I/O Options");
        System.out.println("\t-out    [path]   : Output directory for ROOT file  (def: \$PWD)");
        System.out.println("\t-f               : Overwrite existing files without asking first");
        System.out.println("\t-tree   [string] : Output tree name for ROOT file  (def: \"t\")");
        System.out.println("\t-d      [path]   : Path to file with options");
        System.out.println("\t-log    [path]   : Log file path for analysis");
        System.out.println("\t-n      [int]    : Maximum # of events to add      (def: -1)");
        System.out.println("\t-nf     [int]    : Maximum # of files to analyze   (def: 10)");
        System.out.println("\t-sp     [int]    : Split output after set # events (def: 0)");
        System.out.println("\t-notify [int]    : Print message after set # of events");
        System.out.println("\n * Additional Options:");
        System.out.println("\t-h/--help      : Show this message");
        System.out.println("\t-m/--maps      : Show available pid maps");
        System.out.println("\t-v/--version   : Show version info");
        System.out.println("\t-k/--kin       : Show available kinematics");

        System.out.println("------------------------------------------------------------------------------");
        return false;
    }

    /**
    * Print out available Lund pid to name, charge, mass maps.
    *
    * @param analysis
    *
    * @return Boolean flag
    */
    protected boolean maps(Analysis analysis) {
        System.out.println(" Available PID Maps:");
        System.out.println(" ---------------------------------------------------------------------------- ");
        String[] header = ["Lund PID","Name","Charge (e)","Mass (GeV)"]; System.out.format("  %-15s%-15s%-15s%-15s\n", (Object[])header);
        System.out.println(" ---------------------------------------------------------------------------- ");
        ArrayList<Integer> sortedKeys = new ArrayList<>(analysis.getConstants().getNMap().keySet()); Collections.sort(sortedKeys);
        for (int key : sortedKeys) {
            String[] row = [Integer.toString(key), analysis.getConstants().getName(key), Integer.toString(analysis.getConstants().getCharge(key)), Double.toString(analysis.getConstants().getMass(key))];
            System.out.format("  %-15s%-15s%-15s%-15s\n", (Object[])row);
        }
        System.out.println("------------------------------------------------------------------------------");
        return false;
    }
    
    /**
    * Parse and evaluate command line arguments (or lack thereof).
    * Returns true if main has enough correctly formatted arguments to run smoothly.
    * You will probably want to customize this for less general analyses.
    * @param args      Command line arguments
    * @param analysis  Analysis for which to parse command line arguments
    * @return Boolean flag
    */
    protected boolean parse(String[] args, Analysis analysis) throws IOException {

        System.out.println("------------------------------------------------------------------------------");
        System.out.println(" CLAS12 Analysis \n");

        // version message
        if (args.length==0 || args[0].equals("-v") || args[0].equals("--v") || args[0].equals("-version") || args[0].equals("--version")) { return this.version(); }

		// options message
		if (args[0].equals("--help") || args[0].equals("--h") || args[0].equals("-help") || args[0].equals("-h")) { return this.options(); }

        // pid maps message
        if (args[0].equals("--m") || args[0].equals("--maps") || args[0].equals("-m") || args[0].equals("-maps")) { return this.maps(analysis); }

        // kinematics message
        if (args[0].equals("--k") || args[0].equals("--kin") || args[0].equals("-k") || args[0].equals("-kin")) { return this.kinematics(analysis); }

        // Create log file for recording arguments
        if (args.contains("-log")) {
            Date now = new Date();
            SimpleDateFormat sdf = new SimpleDateFormat ("yyyy.MM.dd.hh.mm.ss");
            String logName = ".analysis."+sdf.format(now)+".log";
            String allArgs = "# CLAS12 Analysis v. 1.0 : Argument list:\n";
            for (String arg : args) { allArgs+=arg; if (!arg.startsWith("-")) { allArgs+=" \n"; } else { allArgs+=" "; } }
            allArgs += "\n";
            System.out.println(" Log file created at: "+logName);
            try {
                File logFile = new File(".analysis"+sdf.format(now)+".log");
                FileWriter writer = new FileWriter(logName);
                writer.write(allArgs); writer.close();
                //analysis.setLogFile(logName); //TODO: pass to analysis object for more extensive logging
            } catch (IOException e) { e.printStackTrace(); return; }
        }

        // Check for drive file
        if (args[0].equals("-d")) {
            String drivePath = args[1];
            try {
                File driveFile = new File(drivePath);
                Scanner reader = new Scanner(driveFile);
                System.out.println(" Loading options:");
                String options = "";
                while (reader.hasNextLine()) {
                    String opt = reader.nextLine();
                    if (opt.startsWith("#")) { continue; } //NOTE: Ignore comments
                    opt = opt.replaceAll("\n"," ");
                    if (!opt.endsWith(" ")) { opt += " "; }
                    options += opt; System.out.println("\t"+opt);
                }
                reader.close();
                args = options.split(" ");
            }
            catch (IOException e) {
                System.out.println(" ERROR: Could not read drive file at: \n "+drivePath);
                System.out.println(" EXITING...");
                System.out.println("------------------------------------------------------------------------------");
                return false;
            }
        }

        // Check input file/directory for HIPO files
        if (!args[0].endsWith("/") && !args[0].endsWith(".hipo")) { args[0]+="/"; }
        analysis.setInPath(args[0]);
        File folder = new File(args[0]);
        if (!folder.isFile() && !folder.isDirectory()) {
            System.out.println(" ERROR: Could not access input files: "+folder);
            System.out.println(" EXITING...");
            System.out.println("------------------------------------------------------------------------------");
            return false;
        }
        if (folder.isFile() && !(folder.getName().endsWith(".hipo"))) {
            System.out.println(" ERROR: Input file: "+folder+" is not a HIPO file.");
            System.out.println(" EXITING...");
            System.out.println("------------------------------------------------------------------------------");
            return false;
        }
        if (folder.isDirectory()) {
            boolean found_hipo = false;
            File[] listOfFiles = folder.listFiles();
            for (File file : listOfFiles) { if (file.isFile() && file.getName().contains(".hipo")) { found_hipo = true; } }
            if (!found_hipo) {
                System.out.println(" ERROR: Input directory: ");
                System.out.println("        "+folder);
                System.out.println("        contains no HIPO files.");
                System.out.println(" EXITING...");
                System.out.println("------------------------------------------------------------------------------");
                return false;
            }
        }

		// Interpret command line arguments
        boolean valid_opt = false;
        for (int i=0; i<args.length; i++) {
            String arg = args[i];

            // Check for simple options in command line
            switch(arg) {

                // Output directory option
                case "-out":
                    if (args.length<=2) { break; }
                    analysis.setOutPath(args[i+1]);
                    File out = new File(analysis.getOutPath());
                    if (out.isDirectory()) {
                        System.out.println(" WARNING: Setting output file name to Analysis.root");
                        if (analysis.getOutPath().endsWith("/")) { analysis.setOutPath(analysis.getOutPath()+"Analysis.root"); }
                        else { analysis.setOutPath(analysis.getOutPath()+"/Analysis.root"); }
                    }
                    // if (out.isFile()) { // Double coded
                    //     System.out.println(" WARNING: "+analysis.getOutPath());
                    //     System.out.print("          already exists, do you want to overwrite? (y/n): ");
                    //     Scanner scanner = new Scanner(System.in);
                    //     String permission = scanner.nextLine();
                    //     if (!(permission.equals("y") || permission.equals("yes") || permission.equals("Yes") || permission.equals("YES"))) {
                    //         System.out.println(" EXITING...");
                    //         System.out.println("------------------------------------------------------------------------------");
                    //         return false;
                    //     }
                    // }
                    valid_opt = true; break;

                // Decay channel specification
                case '-ch':
                    if (args.length<=2) { break; }
                    try {
                        String[] arr = args[i+1].replace(',',':').split(':'); //TODO: Make sure to change later instances!

                        // Get particle groups delimited with just commas
                        ArrayList<ArrayList<String>> arrnew = (ArrayList<ArrayList<String>>)args[i+1].split(':').collect{ el -> return el.split(',')};
                        int k = 0;
                        ArrayList<ArrayList<Integer>> groups = new ArrayList<ArrayList<Integer>>();
                        for (ArrayList<String> group : arrnew) {
                            ArrayList<Integer> list = new ArrayList<Integer>();
                            for (String pid : group) { list.add(k++); }
                            if (list.size()>1) groups.add(list);
                        }

                        // Get list of pids
                        ArrayList<Integer> pids = new ArrayList<Integer>();
                        for (String entry : arr) { int pid = Integer.parseInt(entry); pids.add(pid); }
                        analysis.setDecayAndGroups(pids,groups); //NOTE: Important to use this method since it sorts groups to match the decay.
                        if (args.contains('-mch'))  { analysis.setCombo(true); } // Use processComboEvents() method if both MC::Lund and REC::Particle channel are specified.
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); }

                // Parent channel specification
                case '-pch':
                    if (args.length<=2) { break; }
                    try {
                        String[] arr = args[i+1].split(':');
                        ArrayList<Integer> pids = new ArrayList<Integer>();
                        for (String entry : arr) { int pid = Integer.parseInt(entry); pids.add(pid); }
                        if (pids.size()==0) { return this.help(); }
                        analysis.setParents(pids);
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); }

                // MC channel specification
                case '-mch':
                    if (args.length<=2) { break; }
                    try {
                        String[] arr = args[i+1].replace(',',':').split(':'); //TODO: Make sure to change later instances!

                        // Get particle groups delimited with just commas
                        ArrayList<ArrayList<String>> arrnew = (ArrayList<ArrayList<String>>)args[i+1].split(':').collect{ el -> return el.split(',')};
                        int k = 0;
                        ArrayList<ArrayList<Integer>> groups = new ArrayList<ArrayList<Integer>>();
                        for (ArrayList<String> group : arrnew) {
                            ArrayList<Integer> list = new ArrayList<Integer>();
                            for (String pid : group) { list.add(k++); }
                            if (list.size()>1) groups.add(list);
                        }

                        // Get list of pids
                        ArrayList<Integer> pids = new ArrayList<Integer>();
                        for (String entry : arr) { int pid = Integer.parseInt(entry); pids.add(pid); }
                        if (pids.size()==0) { return this.help(); }
                        analysis.setMCDecay(pids);//TODO: May need set groups and decay here too...
                        if (!args.contains("-ch")) { analysis.setDecayAndGroups(pids,groups); analysis.setUseMC(true); } //NOTE: just use analysis._decay object for simplicity if just looking in MC::Lund bank.  Also, it's important to use this method so that groups is sorted correctly.
                        if (args.contains('-ch'))  { analysis.setCombo(true); } // Use processComboEvents() method if both MC::Lund and REC::Particle channel are specified.
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); }

                case "-q":   analysis.setRequirePID(false); valid_opt = true; break;
                case "-vtx": analysis.setAddVertices(true); valid_opt = true; break;
                case "-ang": analysis.setAddAngles(true); valid_opt = true; break;
                case "-e":   analysis.setRequireE(false); valid_opt = true; break;
                case "-s":   analysis.setStrict(true); valid_opt = true; break;
                case "-rn":  analysis.setAddRunNum(true); valid_opt = true; break;
                case "-en":  analysis.setAddEvNum(true); valid_opt = true; break;
                case "-ml":  analysis.setAddML(true); valid_opt = true; break;
                case "-lk":  analysis.setLambdaKin(true); valid_opt = true; break;
                case "-ak":  analysis.setAffKin(true); valid_opt = true; break;
                case "-ik":  analysis.setIndivKin(true); valid_opt = true; break;
                case "-qa":  analysis.setQA(true); valid_opt = true; break;
                case "-fc":  analysis.setFC(true); valid_opt = true; break;
                case "-ma":  analysis.setMatch(true); valid_opt = true; break;
                // case "-mc":  analysis.setUseMC(true); valid_opt = true; break;

                // Fiducial cuts option
                case "-set_run":
                    if (args.length>2) { try { analysis.setRunNum(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { System.out.println(" WARNING: No run number supplied.  Defaulting to bank values."); } }
            
                // QA method option
                case "-qm":
                    if (args.length>2) { try { analysis.setQAMethod(args[i+1]); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Fiducial cuts option
                case "-fc":
                    if (args.length>2) { try { analysis.setFCLevel(true,Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { System.out.println(" WARNING: Using loose fiducial cuts.  Level not specified."); } }

                // Momentum corrections option
                case "-momc":
                    if (args.length>2) { try { analysis.setMCVersion(args[i+1],Integer.parseInt(args[i+2]),Integer.parseInt(args[i+3])); valid_opt = true; break; }
                    catch (Exception exception) { System.out.println(" WARNING: Using Momentum Correction version for dataset: Fall2018, pass: 1."); } }

                // Momentum corrections option
                case "-mcsmear":
                    if (args.length>2) { try {
                        analysis.setMCSmearing(true);
                        analysis.loadMCSmearingJSON(args[i+1]);
                        analysis.setMCSmearing(Double.parseDouble(args[i+2]),Double.parseDouble(args[i+3]),Double.parseDouble(args[i+4]));
                        analysis.setMCSmearingUseMu(Boolean.parseBoolean(args[i+5]));
                        valid_opt = true; break;
                    } catch (Exception exception) { return this.help(); } }

                // PID tag option
                case "-tag":
                    if (args.length>2) { try {
                        if (!args[i+1].contains(':')) { analysis.setTag(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                        String[] arr = args[i+1].split(':');
                        ArrayList<Integer> pids = new ArrayList<Integer>();
                        for (String entry : arr) { int pid = Integer.parseInt(entry); pids.add(pid); }
                        if (pids.size()==0) { return this.help(); }
                        analysis.setTag(pids);
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); } }

                // Exclusive tag option
                case "-ex":
                    if (args.length>2) { try {
                        analysis.setExclusive(true);
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); } }

                // PID tag option
                case "-filter":
                    if (args.length>2) { try {
                        ArrayList<Integer> pids    = new ArrayList<Integer>();
                        ArrayList<Integer> filters = new ArrayList<Integer>();
                        String[] arr = (String[])(arg.contains(':') ? args[i+1].split(':') : [args[i+1]]);
                        for (String el : arr) {
                            pids.add(Integer.parseInt(el.split(',')[0]));
                            filters.add(Integer.parseInt(el.split(',')[1]));
                        }
                        if (pids.size()==0) { return this.help(); }
                        analysis.setPidFilter(pids,filters);
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); } }

                // ML Client host and port option
                case "-mlclient":
                    if (args.length>2) { try {
                        String host = "localhost";
                        int port = 5000;
                        if (args[i+1].contains(',')) {
                            host = args[i+1].split(',')[0];
                            port = Integer.parseInt(args[i+1].split(',')[1]);
                        } else {
                            host = args[i+1];
                        }
                        analysis.setMLClient(host,port);
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); } }

                // ML Client input bank names option
                case "-mlclient_banks":
                    if (args.length>2) { try {
                        ArrayList<String> banks = new ArrayList<String>();
                        if (args[i+1].contains(',')) {
                            String[] arr = args[i+1].split(',');
                            for (String entry : arr) { banks.add(entry); }
                        } else { banks.add(args[i+1]); }
                        if (banks.size()==0) { return this.help(); }
                        analysis.setMLClientInputBanks(banks);
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); } }

                // ML Client task option
                case "-mlclient_nscores":
                    if (args.length>2) { try {
                        Integer nScores = Integer.parseInt(args[i+1]);
                        analysis.setMLClientNScores(nScores);
                        valid_opt = true; break;
                    }
                    catch (Exception exception) { return this.help(); } }

                // Number of files option
                case "-nf":
                    if (args.length>2) { try { analysis.setNFiles(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Number of events (added) option
                case "-n":
                    if (args.length>2) { try { analysis.setNEvents(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Split number option
                case "-sp":
                    if (args.length>2) { try { analysis.setSplit(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Beam Energy option
                case "-be":
                    if (args.length>2) { try { analysis.setBeamE(Float.parseFloat(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Target Mass option
                case "-tm":
                    if (args.length>2) { try { analysis.setTargetM(Float.parseFloat(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Target Mass option
                case "-tpid":
                    if (args.length>2) { try { analysis.setTargetPID(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Target spin vector option
                case "-tspin":
                    if (args.length>2) {
                        try {
                            double px = Float.parseFloat(args[i+1]);
                            double py = Float.parseFloat(args[i+2]);
                            double pz = Float.parseFloat(args[i+3]);
                            LorentzVector lv_s = new LorentzVector();
                            lv_s.setPxPyPzM(px,py,pz,(double)0.0);
                            analysis.setTargetSpinLV(lv_s);
                            valid_opt = true; break;
                        }
                    catch (Exception exception) { return this.help(); } }

                // Target spin sign option
                case "-tspin_sign":
                    if (args.length>2) { try { analysis.setTSpinSign(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // xF cut option
                case "-xF":
                    if (args.length>2) { try { analysis.setMinxF(Float.parseFloat(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // max mass cut option
                case "-m":
                    if (args.length>2) { try { analysis.setMaxMass(Float.parseFloat(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // notify option
                case "-notify":
                    if (args.length>2) { try { analysis.setNotify(Integer.parseInt(args[i+1])); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                // Tree name option
                case "-tree": 
                    if (args.length>2) { try { analysis.setTreeName(args[i+1]); valid_opt = true; break; }
                    catch (Exception exception) { return this.help(); } }

                default:
                    break;

            } //switch(arg)

            // generic kinematic cut option
            for (String var : analysis.getKinematics().keySet()) {
                if ( arg.startsWith('-'+var+'>') ) {
                    try {
                        System.out.println(' Added cut: '+arg.substring(1));
                        String newstring = new String(arg.substring(2+var.length())); //NOTE: Not really sure why this is necessary for helicity cuts...but it works.
                        Cut cut = (double x) -> {if(x > Float.parseFloat(newstring) ) {return true;} else {return false;}};
                        analysis.addCut(var,cut); valid_opt = true; }
                    catch (Exception exception) { return this.help(); }
                }
                if (arg.startsWith('-'+var+'<') ) {
                    try {
                        System.out.println(' Added cut: '+arg.substring(1));
                        String newstring = new String(arg.substring(2+var.length())); //NOTE: Not really sure why this is necessary for helicity cuts...but it works.
                        Cut cut = (double x) -> {if(x < Float.parseFloat(newstring) ) {return true;} else {return false;}};
                        analysis.addCut(var,cut); valid_opt = true; }
                    catch (Exception exception) { return this.help(); }
                }
            } // for(String var : analysis.getKinematics()...)

        } // loop over args

        // Help when given invalid options
        if (!valid_opt && args.length>1) {
            System.out.println("Invalid option: Use the --help option to see proper usage.");
            System.out.println("------------------------------------------------------------------------------");
            return false;
        }

        // Check if $PWD/Analysis.root file already exists
        File out = new File(analysis.getOutPath());
        if (out.isFile() && !args.contains("-f")) {
            System.out.println(" WARNING: "+analysis.getOutPath());
            System.out.print("          already exists, do you want to overwrite? (y/n): ");
            Scanner scanner = new Scanner(System.in);
            String permission = scanner.nextLine();
            if (!(permission.equals("y") || permission.equals("yes") || permission.equals("Yes") || permission.equals("YES"))) {
                System.out.println(" EXITING...");
                System.out.println("------------------------------------------------------------------------------");
                return false;
            }
        }

        return true;
    } // parseArgs

} // class
