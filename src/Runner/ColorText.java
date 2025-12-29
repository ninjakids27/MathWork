package Runner;

/**
 * This class provides ANSI escape codes for coloring and formatting text in the terminal.
 * It includes constants for various text colors, background colors, and text styles,
 * as well as methods for applying common formatting styles to strings.
 */

public class ColorText {
    // basically have custom colored text in the terminal to have an aesthetically pleasing environment
    // RESET AND TEXT FORMAT
    public static final String CLEAR_SCREEN       = "\u001B[2J";
    public static final String RESET              = "\u001B[0m";
    public static final String BOLD               = "\u001b[1m";
    public static final String UNDERLINE          = "\u001B[4m";
    public static final String ITALIC             = "\u001b[3m";
    public static final String BLINK              = "\u001B[5m";

    // CURSOR CONTROL
    // TDL
    // TEXT COLORS
    public static final String BLACK              = "\u001B[30m";
    public static final String RED                = "\u001B[31m";
    public static final String GREEN              = "\u001B[32m";
    public static final String YELLOW             = "\u001B[33m";
    public static final String BLUE               = "\u001B[34m";
    public static final String PURPLE             = "\u001B[35m";
    public static final String CYAN               = "\u001B[36m";
    public static final String WHITE              = "\u001B[37m";

    // TEXT HIGHLIGHTS
    public static final String BLACK_BACKGROUND   = "\u001B[40m";
    public static final String RED_BACKGROUND     = "\u001B[41m";
    public static final String GREEN_BACKGROUND   = "\u001B[42m";
    public static final String YELLOW_BACKGROUND  = "\u001B[43m";
    public static final String BLUE_BACKGROUND    = "\u001B[44m";
    public static final String PURPLE_BACKGROUND  = "\u001B[45m";
    public static final String CYAN_BACKGROUND    = "\u001B[46m";
    public static final String WHITE_BACKGROUND   = "\u001B[47m";
    
    // METHODS FOR REUSED FORMATTING
    public static String errorFormat(String a){
        return UNDERLINE+BOLD+RED+a+RESET;
    }
    public static String returnFormat(String a){
        return UNDERLINE+BOLD+GREEN+a+RESET;
    }
    public static String dataFormat(String a){
        return UNDERLINE+BOLD+CYAN+a+RESET;
    }
    
    public static void playText(String text, double delay){
        String[] textboxes = text.split(" ");
        for(String S: textboxes){
            System.out.print(S+" ");
            try {
                Thread.sleep((long)delay*1000);
            } catch (Exception e) {
                
            }
        }
        System.out.println();
    }

    public static void clearScreen(){
        System.out.print(CLEAR_SCREEN+RESET);
    }
}
