// Author: Bernard Hollands


using System;
using System.IO;
namespace Perceptron
{
    class Vector{
        private double [] data = null;
        public int Size{get{return data.Length;}}
        public double this[int i]{
            get{return data[i];}
            set { data[i] = value; }
        }

        public Vector(){
            int i = 3;
            data = new double [3]; //default size 3
            for(i = 0; i < 3;i++){
                data[i] = 0;
            }
        }
        public Vector(int size){
            if(size <= 0){
                throw new Exception("Cannot Make Vector Less Than Size 1");
            }
            data = new double[size];
            for(int i = 0; i < size;i++){
                data[i] = 0;
            }
        }
        public Vector(double[] data){
            this.data = new double [data.Length];
            for(int i =0; i < data.Length;i++){
                this.data[i] = data[i];
            }
        }


        public static Vector RandomVector(int length){
            if(length < 1){
                throw new Exception("Cannot Make Vector Less Than Size 1");
            }
            Random rnd = new Random();
            double [] tmp = new double[length];
            for(int i = 0; i < length; i++){
                tmp[i] = rnd.NextDouble();
            }
            Vector rand = new Vector(tmp);
            return rand;
        }

        public static Vector ZeroVector(int length){
            if(length < 1){
                throw new Exception("Cannot Make Vector Less Than Size 1");
            }
            double [] tmp = new double[length];
            for(int i = 0; i < length; i++){
                tmp[i] = 0.0;
            }
            Vector zero = new Vector(tmp);
            return zero;
        }
            
        //dot product
        public static double operator*(Vector left,Vector right){
            if(left.Size != right.Size){
                throw new Exception("Cannot Dot Product Vectors of size "+ left.Size +" and "+ right.Size);
            }
            double tmp = 0;
            for(int i = 0; i < left.Size; i++){
                tmp+= left.data[i]*right.data[i];
            }
            return tmp;
        }
        public static Vector operator*(double left,Vector right){
            Vector tmp = new Vector(right.Size);
            for(int i = 0; i < right.Size; i++){
                tmp.data[i] = right.data[i]*left;
                //tmp.OverWrite(i,);
            }
            return tmp;
        }
        public static Vector operator+(Vector left, Vector right){
            if(left.Size != right.Size){
                throw new Exception("Cannot add Vectors of size "+ left.Size +" and "+ right.Size);
            }
            Vector tmp = new Vector(left.Size);
            for(int i = 0;i < tmp.Size;i++){
                tmp.data[i] = left.data[i] + right.data[i];
            }
            return tmp;
        } 

        /*
        Method to take data from CSV file and put into 2D array of doubles

        INPUT:                           | OUTPUT:  
        toAdd = the vector to append on  | total = new vector with the new vector on the end
                the end of current data  | e.g {1,2,3}.append({4,5,6}) = {1,2,3,4,5,6}
        */
        public Vector Append(Vector toAdd){ //puts one vector onto the end of another
            Vector total = new Vector(data.Length+toAdd.Size); //new vector is length of both combined
            for(int i = 0; i < data.Length; i++){ //for the length of the current vector
                total[i] = data[i]; //simply copy the data in
            }
            for(int i = data.Length; i < total.Size;i++){ //start with i at the the end point of last loop
                total[i] = toAdd[i-data.Length]; //add the input vectors to the total output vector
            }
            return total; //output vector
        }

        /*
        Overloads to string class to print vector
        */
        public override string ToString(){
            string tmp =""; //set temp string
            for(int i = 0;i<Size;i++){ //for the size of the vector
                tmp+=data[i].ToString(); //for each element turn into string and add to tmp
                if(i != Size-1){ //if not at end position
                    tmp+=","; //add a comman
                }  //This is avoid a comma after the final number on printing
            }
            return tmp; //return the string
        }
        public double Get(int index){
            if(index > Size || index < 0){
                throw new Exception("Cannot get value form index that does not exist");
            }
            return data[index];
        }

        public void Set(int index,double value){
            if(index > Size || index < 0){
                throw new Exception("Cannot not overwrite value that does not exist");
            }
            data[index] = value;
        }
       
    }

    /*
    Static Class to do major data handling such as reading for CSV and get methods to easily 
    get 2d array columns and rows and makes perceptron class cleaner.
    */
    static class DataHandeler{
        /*
        Method to take data from CSV file and put into 2D array of doubles

        INPUT:                    | OUTPUT: 
        path = Disk path to file  | dataArray = 2D array with number from CSV without header
        */
        public static double[,] ReadCSV(string path){
            if(!path.Contains(".csv")) //check if filename/path is actually to as csv
                throw new Exception("File "+ path +" is not a csv");
            int[] csvInfo = GetCSVSize(path); //get the length,width and header size of csv file
            double[,] dataArray =  new double[csvInfo[1],csvInfo[2]]; //initalise dataArray with length (csvInfo[1]) and width (csvInfo[2])
            int headerSize = csvInfo[0]; //set the header size
            int row = 0; //set row counter
            StreamReader rdr = new StreamReader(path);
            while(rdr.Peek() >-1){ //while there is a line to read
                if(row < headerSize){ //check if at the header
                    string tmp = rdr.ReadLine(); //read line but do nothing with it
                }else{ //if out of header
                    string tmp = rdr.ReadLine(); //read line into tmp
                    string[] sarr = tmp.Split(','); //split it by the comma
                    for(int i = 0; i < sarr.Length;i++){ //for the length of the split array
                        double.TryParse(sarr[i],out dataArray[row - headerSize,i]); //Add as double to array row "row - headersize"
                    }
                } 
                row++; //increase row after loop
            }
            rdr.Close();//close the file
            return dataArray;//return the data
        }

        /*
        Method to find the size of data in the CSV file including length,width and header size
        Will throw an exception if data in CSV has an error e.g Inconsistant number of colummns per row 
        
        INPUT:                    | OUTPUT: 
        path = Disk path to file  | csvInfo = integer array where each element stores a different file parameter
                                  | csvInfo[0] = Size of CSV file Header
                                  | csvInfo[1] = Number of rows in CSV file
                                  | csvInfo[2] = Number of columns in CSV file */
        public static int[] GetCSVSize(string path){
            StreamReader rdr = new StreamReader(path);
            int length = 0; //set the length of file to 0
            int refWidth = 1; //a valid csv will always have atleast one column
            int headerSize = 0; //set header size to zero
            while(rdr.Peek()>-1){ //while there is a line to be read
                string tmp = rdr.ReadLine(); //read in line from file
                string[] sarr = tmp.Split(','); //split it back on the commas
                if(double.TryParse(sarr[0],out double result)){ //check header by tying to turn into double(All csv header I have ever seen has been letters)
                    length++;
                    if(length == 1) {refWidth = sarr.Length;}//if on first iteration set previous width to current width
                    CheckCellsHaveValues(sarr); //checks that every cell is a number and is present
                    CheckCSVWidth(sarr,refWidth);//checks that the incoming data have consistant widths
                }else{ // if cant be converted 
                    headerSize++; //increase header size
                }
            }
            rdr.Close();
            
            //Create array with all the info collected about the size of data and pass it out 
            int[] csvInfo = {headerSize,length,refWidth};
            return csvInfo;
        }

        /*
        Method to check if each element/cell of the string array(from the CSV) can be converted into a double  
        
        INPUT:                   | OUTPUT: 
        arr = row from csv       | Null
        */
        private static void CheckCellsHaveValues(string[] arr){
            foreach(string str in arr){ //for each stringin the input array
                if(!double.TryParse(str, out double result)){ //check if it can be converted to a double
                    // if fails throw exception
                    throw new Exception("Not all CSV values can be converted to double value, Check source Data");
                }
            }
        }

        /*
        Method to check if each element/cell of the string array(from the CSV) can be converted into a double  
        
        INPUT:                                          | OUTPUT: 
        sarr = row from csv                             | Null
        refWidth = reference width to compare to sarr   |
        */
        private static void CheckCSVWidth(string[] sarr, int refWidth){
            if(sarr.Length != refWidth){//on nex loops if the current width is not same as previos or somthing is missing
                throw new Exception("The CSV columns have changing widths, Check source Data");
            }
        }

        /*
        Method that returns a single column from a given 2d Array
        
        INPUT:                                    | OUTPUT: 
        arr = 2D Array to find column from        | column = 1D array of doubles with values of specified column
        columnIndex = Index of colulm to return   |
        */
        public static double[] GetArrayColumn(double[,] arr, int columnIndex){
            int numRows = arr.GetLength(0);
            double[] column = new double [numRows];
            for(int i = 0; i < numRows; i++){ //for the number of rows 
                column[i] = arr[i,columnIndex];
            }
            return column;
        }

        /*
        Method that gets and returns elements from a row of a 2d Array
        
        INPUT:                                    | OUTPUT: 
        arr = 2D Array                            | output = 1D array with values between start and stop from row
        row = Index of row to get elements from   | e.g row from arr = [1,2,3,4,5]
        start = index of first element to return  | start = 1
        stop =  index of final element to return  | stop = 3
                                                  | output = [2,3,4]*/
        public static double[] Get2DArrayRowSubSection(double[,] arr, int row, int start,int stop){
            int numRows = arr.GetLength(1);
            double[] output = new double [stop];
            for(int i = start; i <= stop; i++){ //for the number of rows 
                output[i - start] = arr[row,i];
            }
            return output;
        }
    
    }
    class Perceptron{
        private Vector weights = null; //vector to hold the weight values
        private Vector[] inputs = null; //array to hold all the input vectors (essentially a matrix)
        private Vector targets = null; // vector to hold the target values
        private int maxIterations; //int to hold the max number iterations can be completed before breaking
        private bool trained = false; //bool to keep track of if perceptorn has been trained
        private int MaxIterations{
            get{return maxIterations;}
            set{
                if(value <= 0){
                    throw new Exception("Max Iterations must be greater than 0");
                }
                maxIterations = value;
            }
        }
        private double learningRate;
        private double LearningRate{
            get{return learningRate;}
            set{
                if(value <= 0){
                    throw new Exception("Learning Rate must be greater than 0");
                }
                learningRate = value;
            }
        }
        /*
        Provide 4 options for initalising perceptron 
        */

        public Perceptron(){
            LearningRate = 0.1;
            this.MaxIterations = 10000000;
        }
        public Perceptron(double learningRate){
            LearningRate = learningRate;
            MaxIterations = 10000000;

        }

        public Perceptron(int maxIterations){
            LearningRate = 0.1;
            MaxIterations = maxIterations;

        }
        public Perceptron(double learningRate,int maxIterations){
            LearningRate = learningRate;
            MaxIterations = maxIterations;

        }

        /*
        Reads data from the CSV into the Perceptron enviroment and gets 
        the inputs vectors, targets and generates the weights
        
        INPUT:                                                          | OUTPUT: 
        filename = The Disk path to the CSV file                        | null
        (only file name is needed if CSV in in same directory as code)
        */
        public void ReadData(string filename){
            double[,]data = DataHandeler.ReadCSV(filename); //read the daat in from the csv skipping 1 line
            GetTargetValues(data,3); //Gets the target values(yi) for the perceptron (the targets are in the 3rd column)
            GetInputVectors(data,1,2); //gets the input values from the data (Input numbers are in column 1 and 2)
            AddBiasToInputs(); //Once inputs are formatted and in array add bias
            GenerateWeights();//Once data is read in we know how big to make weights vector
        }
        
        /*
        Gets the targets from a specifid column in a 2D array of doubles 
        and stores the entire column in a vector
        
        INPUT:                                                       | OUTPUT: 
        data = 2d Array with CSV data                                | null
        column = the index of the column that contraines the targets |
        */
        private void GetTargetValues(double[,] data, int column){
            double[] targetsArr = DataHandeler.GetArrayColumn(data,column); // get the specifed column using the data handler in array
            targets = new Vector(targetsArr); //make targets a vector for consistancy
        }

        /*
        Makes the input vectors and puts them into an array (Matrix),
        Works on assumnption that in CSV file all columns
        with input values will be beside each other
        
        INPUT:                                               | OUTPUT: 
        data = 2D Array with CSV data                        | null
        startColumn = index of first column of input values  |
        startEnd = index of last column of input values      |
        */
        private void GetInputVectors(double[,] data,int startColumn, int endColumn){
            int length = data.GetLength(0); //length of data
            int width = (endColumn - startColumn)+1; //width of columns we want to extract
            inputs = new Vector[length]; //initalise input array to hold vectors

            for(int i = 0; i < length; i++){//for the length of the data
                double[]tmp = new double[width]; //make tmp array to store the data
                tmp = DataHandeler.Get2DArrayRowSubSection(data,i,startColumn,endColumn); //get the row from data between the specified indexes, 
                Vector row = new Vector(tmp); //make those values a vector
                this.inputs[i] = row; //store that vector as an element in the input vector array
            }
        }

        /*
        Adds bias to the input vector, 
        
        INPUT: | OUTPUT: 
        null   | null
        */
        private void AddBiasToInputs(){
            double[] bias = {1}; //make array with bias value 1
            Vector biasVector = new Vector(bias); //make array into vector
            for(int i = 0; i < inputs.Length; i++){ //for the amaount of input vectors 
                inputs[i] = biasVector.Append(inputs[i]); //append the the inputs onto the bias vector and put back in same position
            }
        }

        /*
        Creates the vector of weights same width as input vectors, 
        
        INPUT: | OUTPUT: 
        null   | null
        */
      private void GenerateWeights(){
            weights = Vector.ZeroVector(inputs[0].Size); //make a weights vector of same as size as the inputs
        }

        /*
        Trains the weights of the perceptron using the Gradient Descent algorithum, 
        
        INPUT: | OUTPUT: 
        null   | null
        */
        public void TrainData(){
            if(inputs == null || targets == null){ //Check that inputs and targets exist
                throw new Exception("There was no data provided for Perceptron to train on!");
            }
            int error; //initalise errors
            int itr = 0; //set iterations
            do{ 
                error = 0; //Step 1 set error to zero
                for(int i = 0; i < inputs.Length; i++){ //Step 2
                    double yhat = StepFunction(inputs[i]*weights); //get estimation of input vector dot weights
                    if(yhat != targets[i]){ //if yhat is not correct 
                        weights = weights + learningRate *(targets[i] - yhat)*inputs[i]; //update weights
                        error++; //increase error
                    }
                }
                itr++; //increase iteration count
                if(itr == maxIterations){ //check if max iterations have been reached
                    throw new Exception("Perceptron Failed to Converge after "+ maxIterations +" iterations");} //tell user that perceptron did not converge
            }while(error > 0); //if there was an error go back to step 1 otherwise move on
            trained = true; //acknoledge that training was succsessful
            Console.WriteLine("Perceptron Converged after {0} Iterations",itr); //tell user that perceptron converged
        }   
        
        /*
        Simple step function that outputs 1 if x > 0 
        
        INPUT:            | OUTPUT: 
        x = input value   | 1 or 0       
        */
        public double StepFunction(double x){
            if (x <= 0){ //if x less than or = to zero
                return 0;
            }
            return 1; //if x is larger than 0
        }
        /*
        Shows formatted version of what the perceptron has calculated
        
        INPUT: | OUTPUT: 
        null   | null       
        */
        public void Output(){
            if(!trained){ //check if perceptron has actually been trained
                throw new Exception("The perceptron must be trained before values can be output");
            } //output the weights aswwell as the intercept and slope of found line
            Console.WriteLine("*************************");
            Console.WriteLine("Trained Weight Values are:");
            double w0 = weights[0]; //get each weight value
            double w1 = weights[1];
            double w2 = weights[2];
            Console.WriteLine("w0 = {0}",w0); //print each weight value
            Console.WriteLine("w1 = {0}",w1);
            Console.WriteLine("w2 = {0}",w2);
            Console.WriteLine("*************************");
            Console.WriteLine("Found Line has:");
            double intrcept = -(w0)/(w2); //calculate intercept
            Console.WriteLine("Intercept = {0}",intrcept);
            double slope = -(w0/w2)/(w0/w1); //calculuate slope
            Console.WriteLine("Slope = {0}",slope);
            Console.WriteLine("*************************");
        }

        /*
        Calculates output for a single input vector
        
        INPUT:                                            | OUTPUT: 
        inputVec = input vector with values to classify   | result = classification result       
        */
        public int ClassifyPoint(Vector inputVec){ 
            if(!trained){ //check if perceptron is trained
                throw new Exception("The perceptron must be trained before values can be output");
            }
            if(inputVec.Size != 2){ //check vector is right size
                throw new Exception("Input Vector must be of length 2");
            }
            double[] arr = {1}; //make bias array
            Vector bias = new Vector(arr); //make bias vector
            Vector data = bias.Append(inputVec); //add bias to input Vector
            int result = (int)StepFunction(data*weights); //calculate output result
            return result;
        }
    
    }
}
