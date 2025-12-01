import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
public class StatOps {
    
    private static final double epsilon = 1e-9;
    private static final double bigEpsilon = 1e15;
    /**
     * Uses the defined epsilon for tolerance checking and rounds the input to avoid floating point errors
     * @param input
     * @return
     */
    private static double tolerance(double input){
        double temp = (input*bigEpsilon);
        temp = Math.round(temp);
        return (double)temp/bigEpsilon;
    }
    /**
     * Checks if two doubles are equal within a defined tolerance (epsilon).
     * @param a first double
     * @param b second double
     * @return true if the numbers are considered equal, false otherwise
     */
    public static boolean isEqual(double a, double b) {
        return Math.abs(a - b) < epsilon;
    }

    /**
     * Calculates the mean of an array of doubles.
     * @param vector the array of doubles
     * @return the mean value
     */
    public static double mean(double[] vector){
        double temp = sum(vector);
        temp /= vector.length;
        return temp;
    }

    /**
     * Sums an array of doubles and returns the result as a double.
     * @param array the array to sum
     * @return the sum of the array elements as a double
     */
    public static double sum(double[] array){
        double sum = 0;
        for(int i = 0; i < array.length; i++){
            sum += array[i];
            // check for floating point errors 
            sum = tolerance(sum);
        }
        return sum;
    }

    /**
     * Sums an array of doubles and returns the result as an integer.
     * @param array the array to sum
     * @return the sum of the array elements as an integer
     */
    public static int  sumD(double[] array){
        int sum = 0;
        for(int i = 0; i < array.length; i++){
            sum += array[i];
        }
        return sum;
    }

    /**
     * Calculates the weighted mean of a vector given the weights.
     * @param vector the data vector
     * @param weights the weights for each element
     * @param debug a flag for debugging
     * @return the weighted mean
     */
    public static double weighted_mean(double[] vector, double[] weights,boolean debug){
        double A = MatrixOps.dot_product(vector, weights,debug);
        double B = sum(weights);
        double C = A/B;
        return C;
    }


    /**
     * Sorts the given array in ascending order using the bubble sort algorithm.O(n^2) time complexity.
     * @param arr the array to sort
     * @return the sorted array
     */
    public static double[] sort(double[] arr) {
        // Create a new array to avoid modifying the original
        double[] sortedArr = Arrays.copyOf(arr, arr.length);

        int n = sortedArr.length;
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - i - 1; j++) {
                if (sortedArr[j] > sortedArr[j + 1]) {
                    // Swap sortedArr[j] and sortedArr[j+1]
                    double temp = sortedArr[j];
                    sortedArr[j] = sortedArr[j + 1];
                    sortedArr[j + 1] = temp;
                }
            }
        }
        return sortedArr;
    }

    /**
     * Calculates the median of the given data vector.
     * @param vector the data vector
     * @return the median value
     */
    public static double median(double[] vector){
        double[] sort_list = sort(vector);
        if(sort_list.length % 2 == 0){
            int a = sort_list.length/2;
            int b = sort_list.length/2-1;
            double[] input = {sort_list[a],sort_list[b]};
            return mean(input);
        }else{
            int a = sort_list.length/2;
            return sort_list[a];
        }
    }

    
    /**
     * Generates a random integer between the specified low and high bounds (inclusive).
     * @param lowBound the lower bound (inclusive)
     * @param highBound the upper bound (inclusive)
     * @return a random integer between lowBound and highBound
     */
    public static int random(int lowBound, int highBound){
        return (int) Math.floor(Math.random()*(highBound-lowBound+1)+lowBound);
    }
    /**
     * Generates a frequency list based on the provided data vector and frequency counts.
     * @param vector1 the data vector
     * @param vector2 the frequency counts
     * @return the generated frequency list
     */
    public static double[] freqList(double[] vector1, double[] vector2){
        double[] tempVector = new double[sumD(vector2)];
        int idx = 0;
        for(int i = 0; i < vector2.length; i++){
            for(int k = 0; k < vector2[i]; k++){
                tempVector[idx++] = vector1[i];
            }
        }
        return tempVector;
    }
    

    

    
    /**
     * Finds the indexes of the maximum value(s) in the given vector.
     * @param vector the data vector (int or double type array is fine)
     * @return an array of indexes corresponding to the maximum value(s)
     */
    public static int[] findMax(double[] vector){
        double Max = vector[0];
        List<Integer> indexes = new ArrayList<>();
        for(int i = 0; i < vector.length; i++){
            if(Max < vector[i]){
                Max = vector[i];
                indexes.clear();
                indexes.add(i);
            }else if(Max == vector[i]){
                indexes.add(i);
            }
        }
        int[] return_indexes = new int[indexes.size()];
        for(int i = 0; i < return_indexes.length;i++){
            return_indexes[i] = indexes.get(i);
        }
        return return_indexes;

    }


    public static int[] findMax(int[] vector){
        double Max = vector[0];
        List<Integer> indexes = new ArrayList<>();
        for(int i = 0; i < vector.length; i++){
            if(Max < vector[i]){
                Max = vector[i];
                indexes.clear();
                indexes.add(i);
            }else if(Max == vector[i]){
                indexes.add(i);
            }
        }
        int[] return_indexes = new int[indexes.size()];
        for(int i = 0; i < return_indexes.length;i++){
            return_indexes[i] = indexes.get(i);
        }
        return return_indexes;

    }

    
    public static double[] mode(double[] vector){
        // I decided to kill this shit before it got WAYYY to out of hand
        List<Double> item = new ArrayList<>();
        List<Integer> amount = new ArrayList<>();
        // could do a 2d array but I'm honestly too lazy for that shit

        for(int i = 0; i < vector.length; i++){
            if(item.contains(vector[i])){
                amount.set(item.indexOf(vector[i]), amount.get(item.indexOf(vector[i]))+1);
            }else{
                item.add(vector[i]);
                amount.add(1);
            };
        }
        int[] amount_caveman = new int[amount.size()];
        for(int i = 0; i < amount.size(); i++){
            amount_caveman[i] = amount.get(i);
        }
        int[] a = findMax(amount_caveman);
        double[] b = new double[a.length];
        for(int i = 0; i < a.length;i++){
            b[i] = item.get(a[i]);
        }
        return b;
    }

    /**
     * Calculates the mode of the given data vector.
     * @param vector the data vector (int or double type array is fine)
     * @return the mode value as a double array (to account for multiple modes)
     */
    public static double[] mode(int[] vector){
        // I decided to kill this shit before it got WAYYY to out of hand
        List<Integer> item = new ArrayList<>();
        List<Integer> amount = new ArrayList<>();
        // could do a 2d array but I'm honestly too lazy for that shit

        for(int i = 0; i < vector.length; i++){
            if(item.contains(vector[i])){
                amount.set(item.indexOf(vector[i]), amount.get(item.indexOf(vector[i]))+1);
            }else{
                item.add(vector[i]);
                amount.add(1);
            };
        }
        int[] amount_caveman = new int[amount.size()];
        for(int i = 0; i < amount.size(); i++){
            amount_caveman[i] = amount.get(i);
        }
        int[] a = findMax(amount_caveman);
        double[] b = new double[a.length];
        for(int i = 0; i < a.length;i++){
            b[i] = item.get(a[i]);
        }
        return b;
    }

    /**
     * Calculates the mid-range of the given data vector.
     * @param vector the data vector
     * @return the mid-range value
     */
    public static double midRange(double[] vector){
        double[] a = sort(vector);
        return (a[0]+a[a.length-1])/2;
    }

    /**
     * Calculates the range of the given data vector.
     * @param vector the data vector
     * @return the range value
     */
    public static double range(double[] vector){
        double[] a = sort(vector);
        return (a[a.length-1]-a[0]);
    }

    /**
     * Calculates the first quartile (Q1) of the given data vector.
     * @param vector the data vector
     * @return the first quartile (Q1) value
     */
    public static double getQ1(double vector[]){
        double[] data = sort(vector);
        if(data.length+1 % 4 == 0){
            return data[(data.length+1)/4];
        }else{
            int n1 = (int) (data.length+1.5)/4;
            int n2 = (int) (data.length+0.5)/4;
            double[] b = {vector[n1],vector[n2]};
            double a = mean(b);
            return a;
        }
    }

    /**
     * Calculates the third quartile (Q3) of the given data vector.
     * @param vector the data vector
     * @return the third quartile (Q3) value
     */
    public static double getQ3(double vector[]){
        double[] data = sort(vector);
        if((3*(data.length+1)) % 4 == 0){
            return data[(3*(data.length+1))/4 -1];
        }else{
            int n1 = (int) ((3*(data.length+1)+0.5)/4);
            int n2 = (int) ((3*(data.length+1)-0.5)/4);
            double[] b = {data[n1-1],data[n2-1]};
            double a = mean(b);
            return a;
        }
        
    }

    /**
     * Calculates and prints the five-number summary (minimum, Q1, median, Q3, maximum) for a frequency distribution.
     * @param vector1 data values 
     * @param vector2 frequency counts
     */
    public static void five_number_summary_of_freq(double[] vector1, double[] vector2){
        List<Double> table = new ArrayList<>();
        for(int i = 0; i < vector1.length; i++){
            for(int k = 0; k < (int)vector2[i]; k++){
                table.add(vector1[i]);
            }
        }
        
        double[] table_double = new double[table.size()];
        for (int i = 0; i < table.size(); i++) {
            table_double[i] = table.get(i);
        }
        double[] D = sort(table_double);
        System.out.println("Min: "+D[0]);
        System.out.println("Q1: "+getQ1(D));
        System.out.println("Median: "+median(D));
        System.out.println("Q3: "+getQ3(D));
        System.out.println("Max: "+D[D.length-1]);
    }
    
    
    private static double rounding(double num, int decimal){
        return Math.round(Math.pow(10, decimal)*num)/Math.pow(10, decimal);
    }
    
    

    /**
     * Detects and prints outliers in the given vector using the 1.5*IQR rule.
     * @param vector the data vector
     * @param isSample whether the data is a sample or population
     * @param decimal the number of decimal places to round to
     */
    public static void outlier_test(double[] vector){
        double q1 = getQ1(vector);
        double q3 = getQ3(vector);
        double IQR =  q3 - q1;
        double low_bound = q1-1.5*IQR;
        double high_bound = q3+1.5*IQR;
        for(int i = 0; i < vector.length; i++){
            if(vector[i] > high_bound){
                System.out.println("High outlier found at "+i+" outlier is "+vector[i]);
            }
            if(vector[i] < low_bound){
                System.out.println("Low outlier found at "+i+" outlier is "+vector[i]);
            }
        }
    }


    public static double stdDeviation(double[] vector){
        double mean = mean(vector);
        vector = MatrixOps.vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = MatrixOps.vectorPow(vector, 2);
        double val = sum(vector);
        return Math.sqrt(val/(vector.length));
    }


    public static double stdDeviation(double[] vector,int decimal){
        double mean = mean(vector);
        vector = MatrixOps.vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = MatrixOps.vectorPow(vector, 2);
        double val = sum(vector);
        return rounding(Math.sqrt(val/(vector.length)),decimal);
    }


    public static double samDeviation(double[] vector){
        double mean = mean(vector);
        vector = MatrixOps.vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = MatrixOps.vectorPow(vector, 2);
        double val = sum(vector);
        return Math.sqrt(val/(vector.length-1));
    }


    public static double samDeviation(double[] vector,int decimal){
        double mean = mean(vector);
        vector = MatrixOps.vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = MatrixOps.vectorPow(vector, 2);
        double val = sum(vector);
        return rounding(Math.sqrt(val/(vector.length-1)), decimal);
    }


    public static double samVariance(double[] vector){
        double mean = mean(vector);
        vector = MatrixOps.vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = MatrixOps.vectorPow(vector, 2);
        double val = sum(vector);
        return val/(vector.length-1);
        
    }


    public static double popVariance(double[] vector){
        double mean = mean(vector);
        vector = MatrixOps.vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = MatrixOps.vectorPow(vector, 2);
        double val = sum(vector);
        return val/(vector.length);
        
    }


    public static double popVariance(double[] vector, int decimal){
        double mean = mean(vector);
        vector = MatrixOps.vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = MatrixOps.vectorPow(vector, 2);
        double val = sum(vector);
        return rounding(val/(vector.length), decimal);
        
    }


    public static double Ztest(double observedVal,double[] vector,boolean isSample){
        double z;
        if(isSample){
           z = (observedVal-mean(vector))/samDeviation(vector); 
        }else{
            z = (observedVal-mean(vector))/stdDeviation(vector);
        }
        return z;
    }

    public static double Ztest(double observedVal,double[] vector){
        double z;
        z = (observedVal-mean(vector))/stdDeviation(vector);
        return z;
    }

    public static double Ztest(double observedVal,double mean, double stdDeviation){
           double z = (observedVal-mean)/stdDeviation;
           return z;
    }

    public static double Ztest(double observedVal,double[] vector,boolean isSample,int decimal){
        double z;
        if(isSample){
           z = (observedVal-mean(vector))/samDeviation(vector); 
        }else{
            z = (observedVal-mean(vector))/stdDeviation(vector);
        }
        return rounding(z, decimal);
    }

    public static double Ztest(double observedVal,double[] vector,int decimal){
        double z;
        z = (observedVal-mean(vector))/stdDeviation(vector);
        return rounding(z, decimal);
    }
    /**
     *  Z-test calculation with rounding.
     * @param observedVal the observed value
     * @param vector the data vector it ignores the mean and std deviation
     * @param mean the population mean
     * @param stdDeviation the population standard deviation
     * @param isSample whether the data is a sample or population
     * @param decimal the number of decimal places to round to
     * @return the Z-test statistic
     */
    public static double Ztest(double observedVal,double mean, double stdDeviation,int decimal){
           double z = (observedVal-mean)/stdDeviation;
           return rounding(z, decimal);
    }

    public static void five_number_summary(double[] vector){
        System.out.println("5 number summary:");
        double[] D = sort(vector);
        System.out.println("Min: "+D[0]);
        System.out.println("Q1: "+getQ1(D));
        System.out.println("Median: "+median(D));
        System.out.println("Q3: "+getQ3(D));
        System.out.println("Max: "+D[D.length-1]);
    }


    public static void centralTedencySummary(double[] vector){
        System.out.println("3 measures of central tedency:");
        System.out.println("Mean: "+mean(vector));
        System.out.println("Median: "+median(vector));
        System.out.print("Mode: ");
        MatrixOps.Print_Vector(mode(vector));
    }


    
    
    public static void OneVarStat(double[] vector){
        System.out.println("Mean: "+mean(vector));
        System.out.println("Median: "+median(vector));
        System.out.println("Mode: "+mode(vector));
        System.out.println("Sum: "+sum(vector));
        System.out.println("Sum^2: "+MatrixOps.vectorPow(vector, 2));
        System.out.println("Sam deviation: "+samDeviation(vector));
        System.out.println("Std Deviation: "+stdDeviation(vector));
        System.out.println("Length: "+vector.length);
        five_number_summary(vector);
        
    }
    public static double mean_of_probability_distribution(double[] x, double[] p){
        if(sum(p) != 1){
            throw new IllegalArgumentException(ColorText.errorFormat("Probabilities do not sum to 1"));
        }
        return weighted_mean(x, p,false);
    }
    
    public static double variance_of_probability_distribution(double[] x, double[] p){
        if(sum(p) != 1){
            throw new IllegalArgumentException(ColorText.errorFormat("Probabilities do not sum to 1"));
        }
        double mean = mean_of_probability_distribution(x, p);
        double[] x_minus_mean = MatrixOps.vectorAdd(Arrays.copyOf(x, x.length), -mean);
        double[] x_minus_mean_squared = MatrixOps.vectorPow(x_minus_mean, 2);
        return weighted_mean(x_minus_mean_squared, p,false);
    }
    
    public static double std_of_probability_distribution(double[] x, double[] p){
        return Math.sqrt(variance_of_probability_distribution(x, p));
    }
    
    public static void centralTedencySummary(double[] vector,int decimal){
        System.out.println("3 measures of central tedency:");
        System.out.println("Mean: "+rounding(mean(vector),decimal));
        System.out.println("Median: "+rounding(median(vector),decimal));
        System.out.print("Mode: ");
        MatrixOps.Print_Vector(MatrixOps.roundingVector(mode(vector),decimal));
    }
    public static double mean_of_probability_distribution(double[] x, double[] p,int decimal){
        if(sum(p) != 1){
            throw new IllegalArgumentException(ColorText.errorFormat("Probabilities do not sum to 1"));
        }
        return rounding(weighted_mean(x, p,false), decimal);
    }

    /**
     * Calculates the variance of a probability distribution with rounding.
     * @param x values of the random variable
     * @param p probabilities associated with each x value
     * @param decimal number of decimal places to round to
     * @return variance of the probability distribution
     */
    public static double variance_of_probability_distribution(double[] x, double[] p,int decimal){
        if(sum(p) != 1){
            throw new IllegalArgumentException(ColorText.errorFormat("Probabilities do not sum to 1"));
        }
        double mean = mean_of_probability_distribution(x, p);
        double[] x_minus_mean = MatrixOps.vectorAdd(Arrays.copyOf(x, x.length), -mean);
        double[] x_minus_mean_squared = MatrixOps.vectorPow(x_minus_mean, 2);
        return rounding(weighted_mean(x_minus_mean_squared, p,false), decimal);
    }
    /**
     * Calculates the standard deviation of a probability distribution with rounding.
     * @param x values of the random variable
     * @param p probabilities associated with each x value
     * @param decimal number of decimal places to round to
     * @return standard deviation of the probability distribution
     */
    public static double std_of_probability_distribution(double[] x, double[] p,int decimal){
        return rounding(Math.sqrt(variance_of_probability_distribution(x, p)), decimal);
    }


    public static void probability_distribution(double[] x, double[] p){
        System.out.println("Probability Distribution:");
        System.out.println("Mean: "+mean_of_probability_distribution(x, p));
        System.out.println("Variance: "+variance_of_probability_distribution(x, p));
        System.out.println("Standard Deviation: "+std_of_probability_distribution(x, p));
    }

    /**
     * Prints the probability distribution summary including mean, variance, and standard deviation with rounding.
     * @param x values of the random variable
     * @param p probabilities associated with each x value
     * @param decimal number of decimal places to round to
     */
    public static void probability_distribution(double[] x, double[] p,int decimal){
        System.out.println("Probability Distribution:");
        System.out.println("Mean: "+mean_of_probability_distribution(x, p, decimal));
        System.out.println("Variance: "+variance_of_probability_distribution(x, p, decimal));
        System.out.println("Standard Deviation: "+std_of_probability_distribution(x, p, decimal));
    }


    /**
     * Calculates the binomial probability density function (PDF) for given parameters.
     * @param n number of trials (non-negative integer)
     * @param p probability for success on each trial (0 <= p <= 1)
     * @param x number of successes (0 <= x <= n) it can also be in a list for atleast problems
     * @param individual if true then return the individual probabilities for each x in the list, otherwise return the total probability
     * @return the probability of exactly x successes in n trials
     */
    public static double binomialPDF(int n, double p, int x){
        double q = 1-p;
        double result = numberTheory.nCr(n, x)*Math.pow(p,x)*Math.pow(q,n-x);
        return tolerance(result);
    }

    public static double binomialPDF(int n, double p, int x,int decimal){
        return rounding(binomialPDF(n, p, x), decimal);
    }
    public static double binomialPDF(int n, double p, int[] x){
        double total = 0;
        for(int i = 0; i < x.length; i++){
            total += binomialPDF(n, p, x[i]);
        }
        return tolerance(total);
    }
    public static double binomialPDF(int n, double p, int[] x,int decimal){
        return rounding(binomialPDF(n, p, x), decimal);
    }

    public static double[] binomialPDF(int n, double p, int[] x,boolean individual){
        if(individual){
            double[] results = new double[x.length];
            for(int i = 0; i < x.length; i++){
                results[i] = tolerance(binomialPDF(n, p, x[i]));
            }
            return results;
        }else{
            return new double[]{tolerance(binomialPDF(n, p, x))};
        }
    }
    
    public static double[] binomialPDF(int n, double p, int[] x,boolean individual,int decimal){
        if(individual){
            double[] results = new double[x.length];
            for(int i = 0; i < x.length; i++){
                results[i] = rounding(tolerance(binomialPDF(n, p, x[i])), decimal);
            }
            return results;
        }else{
            return new double[]{rounding(tolerance(binomialPDF(n, p, x)), decimal)};
        }
    }
    /**
     * Calculates the binomial cumulative distribution function (CDF) for given parameters.
     * @param n number of trials (non-negative integer)
     * @param p probability for success on each trial (0 <= p <= 1)
     * @param x number of successes (0 <= x <= n)
     * @return the cumulative probability of up to x successes in n trials
     */
    public static double binomialCDF(int n, double p, int x){
        double total = 0;
        for(int i = 0; i <= x; i++){
            total += tolerance(binomialPDF(n, p, i));
        }
        return total;
    }
    /**
     * Calculates the binomial cumulative distribution function (CDF) for given parameters.
     * @param n number of trials (non-negative integer)
     * @param p probability for success on each trial (0 <= p <= 1)
     * @param x number of successes (0 <= x <= n)
     * @param decimal number of decimal places to round to
     * @return the cumulative probability of up to x successes in n trials
     */
    public static double binomialCDF(int n, double p, int x,int decimal){
        return rounding(binomialCDF(n, p, x), decimal);
    }
    /**
     * Calculates the binomial cumulative distribution function (CDF) for given parameters with options for greater than and at most.
     * @param n number of trials (non-negative integer)
     * @param p probability for success on each trial (0 <= p <= 1)
     * @param x number of successes (0 <= x <= n)
     * @param greaterThan if true, calculates the probability of more than x successes; if false, calculates the probability of up to x successes
     * @param atmost if true, calculates the probability of at most x successes; if false, calculates the probability of exactly x successes
     * @return the calculated probability based on the specified conditions
     */
    public static double binomialCDF(int n, double p, int x,boolean greaterThan){
        if(greaterThan){
            return tolerance(1-binomialCDF(n, p, x));
        }else{
            return tolerance(binomialCDF(n, p, x));
        }
    }
    /**
     * Calculates the binomial cumulative distribution function (CDF) for given parameters with options for greater than and at most, with rounding.
     * @param n number of trials (non-negative integer)
     * @param p probability for success on each trial (0 <= p <= 1)
     * @param x number of successes (0 <= x <= n)
     * @param greaterThan if true, calculates the probability of more than x successes; if false, calculates the probability of up to x successes
     * @param atmost if true, calculates the probability of at most x successes; if false, calculates the probability of exactly x successes
     * @param decimal number of decimal places to round to
     * @return the calculated probability based on the specified conditions, rounded to the specified number of decimal places
     */
    public static double binomialCDF(int n, double p, int x,boolean greaterThan,boolean atmost){
        if(greaterThan){
            return tolerance(1-binomialCDF(n, p, x));
        }
        return tolerance(binomialCDF(n, p, x));
    }
    /**
     * Calculates the binomial cumulative distribution function (CDF) for given parameters with an option for greater than, with rounding.
     * @param n number of trials (non-negative integer)
     * @param p probability for success on each trial (0 <= p <= 1)
     * @param x number of successes (0 <= x <= n)
     * @param greaterThan if true, calculates the probability of more than x successes; if false, calculates the probability of up to x successes
     * @param decimal number of decimal places to round to
     * @return the calculated probability based on the specified condition, rounded to the specified number of decimal places
     */
    public static double binomialCDF(int n, double p, int x,boolean greaterThan,boolean atmost,int decimal){
        return rounding(binomialCDF(n, p, x, greaterThan, atmost), decimal);
    }

    /**
     * Calculates the binomial cumulative distribution function (CDF) for given parameters with an option for greater than, with rounding.
     * @param n number of trials (non-negative integer)
     * @param p probability for success on each trial (0 <= p <= 1)
     * @param x number of successes (0 <= x <= n)
     * @param greaterThan if true, calculates the probability of more than x successes; if false, calculates the probability of up to x successes
     * @param decimal number of decimal places to round to
     * @return the calculated probability based on the specified condition, rounded to the specified number of decimal places
     */
    public static double binomialCDF(int n, double p, int x,boolean greaterThan,int decimal){
        return rounding(binomialCDF(n, p, x, greaterThan), decimal);
    }

    

    public static double normalDistrubution(double mean, double stdDeviation, double x){
        return (1)/(stdDeviation*Math.sqrt(2*Math.PI)*Math.exp(.5*Math.pow((x-mean)/stdDeviation,2)));
    }

    public static double sampleNormalDistrubution(double mean, double stdDeviation){
        return normalDistrubution(mean, stdDeviation, Math.random());
    }
    /**
     * does as described reference notes.pdf for indepth stuff that I don't want to explain
     * @return random z number on the normal distribution
     */
    public static double box_Muller_Transfrom(){
        double u1 = Math.random();
        double u2 = Math.random();
        double s = Math.pow(u1, 2) + Math.pow(u2, 2);
        while(s == 0 || s >= 1){
            u1 = Math.random();
            u2 = Math.random();
            s = Math.pow(u1, 2) + Math.pow(u2, 2);
        }
        return u1 * Math.sqrt((-2*Math.log(s))/s);
    }

    public static double box_Muller_Transfrom(double test1, double test2){
        double u1 = test1;
        double u2 = test2;
        double s = Math.pow(u1, 2) + Math.pow(u2, 2);
        return u1 * Math.sqrt((-2*Math.log(s))/s);
    }
}

