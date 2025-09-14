import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
public class StatOps {
    public static double mean(double[] vector){
        double temp = sum(vector);
        temp /= vector.length;
        return temp;
    }

    
    public static double sum(double[] array){
        double sum = 0;
        for(int i = 0; i < array.length; i++){
            sum += array[i];
        }
        return sum;
    }


    public static int  sumD(double[] array){
        int sum = 0;
        for(int i = 0; i < array.length; i++){
            sum += array[i];
        }
        return sum;
    }


    public static double weighted_mean(double[] vector, double[] weights,boolean debug){
        double A = MatrixOps.dot_product(vector, weights,debug);
        double B = sum(weights);
        double C = A/B;
        return C;
    }


    // Standard bubblesort
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


    public static double[] vectorAdd(double[] vector, double factor){
        for(int i = 0; i<vector.length; i++){
            vector[i] += factor;
        }
        return vector;
    }


    public static double[] vectorMult(double[] vector, double factor){
        for(int i = 0; i<vector.length; i++){
            vector[i] *= factor;
        }
        return vector;
    }


    public static double[] vectorPow(double[] vector, double factor){
        for(int i = 0; i<vector.length; i++){
            vector[i] = Math.pow(vector[i], factor);
        }
        return vector;
    }


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


    public static double midRange(double[] vector){
        double[] a = sort(vector);
        return (a[0]+a[a.length])/2;
    }


    public static double range(double[] vector){
        double[] a = sort(vector);
        return (a[a.length]-a[0]);
    }


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
    public static double[] roundingVector(double[] num, int decimal){
        double[] roundedVector = new double[num.length];
            for(int i = 0; i < num.length; i++){
                roundedVector[i] = rounding(num[i],decimal);
            }
            return roundedVector;
    }

    
    public static void outlier_test(double[] vector){
        
    }


    public static double stdDeviation(double[] vector){
        double mean = mean(vector);
        vector = vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = vectorPow(vector, 2);
        double val = sum(vector);
        return Math.sqrt(val/(vector.length));
    }


    public static double stdDeviation(double[] vector,int decimal){
        double mean = mean(vector);
        vector = vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = vectorPow(vector, 2);
        double val = sum(vector);
        return rounding(Math.sqrt(val/(vector.length)),decimal);
    }


    public static double samDeviation(double[] vector){
        double mean = mean(vector);
        vector = vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = vectorPow(vector, 2);
        double val = sum(vector);
        return Math.sqrt(val/(vector.length-1));
    }


    public static double samDeviation(double[] vector,int decimal){
        double mean = mean(vector);
        vector = vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = vectorPow(vector, 2);
        double val = sum(vector);
        return rounding(Math.sqrt(val/(vector.length-1)), decimal);
    }


    public static double samVariance(double[] vector){
        double mean = mean(vector);
        vector = vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = vectorPow(vector, 2);
        double val = sum(vector);
        return val/(vector.length-1);
        
    }


    public static double popVariance(double[] vector){
        double mean = mean(vector);
        vector = vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = vectorPow(vector, 2);
        double val = sum(vector);
        return val/(vector.length);
        
    }


    public static double popVariance(double[] vector, int decimal){
        double mean = mean(vector);
        vector = vectorAdd(Arrays.copyOf(vector, vector.length), -mean);
        vector = vectorPow(vector, 2);
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


    public static void centralTedencySummary(double[] vector,int decimal){
        System.out.println("3 measures of central tedency:");
        System.out.println("Mean: "+rounding(mean(vector),decimal));
        System.out.println("Median: "+rounding(median(vector),decimal));
        System.out.print("Mode: ");
        MatrixOps.Print_Vector(roundingVector(mode(vector),decimal));
    }
    

    public static void OneVarStat(double[] vector){
        
    }
}
