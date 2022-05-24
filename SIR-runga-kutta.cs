using System;
using System.IO;
/*
Author: Bernard Hollands

-------------------------------------------------------------------------------------------
SIR equations solved with runge-kutta 4th order equations

Example of Hong Kong Flu spreading in new work as found at:
https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model

For this example starting values are:
Situation Information: 
Population/Suseptipal =  7900000
Ininital Infected = 10
R0 = 1.5
D = 3

Initial Values:
suseptipal = 1 (entire population)
infected = 1.26e-6 (10/7900000 = 10 people scaled)
recovered = 0
---------------------------------------------------------------------------------------------

Program.cs file Contains 4 Classes
 Program - Contains main method and driver code
 SIR_Model - Impliments and has functions for Suseptipal, infected and Recovered Model
 RungeKutta - Impliments 4th order Runge kutta method for solving differential equaion systems
 Vector - Enables creation of vectors as storage like arrays and allows vector operations
*/
namespace SIR_solver
{
    //make delegate for funciton that takes a time and vector a dependants
    public delegate double functions(Vector deps);
    class Program
    {
        static void Main(string[] args)
        {
            try{
                double D = 3; //set D
                double R0 = 1.5; //set R0
                double gamma = 3;
                double beta = 3;
                SIR_Model SIR = new SIR_Model(D,R0); //define a SIR modle

                functions[] diffEquations = new functions[3];// = {SIR.DSdt,SIR.DIdt,SIR.DRdt}; //create vector of functions to solve
                diffEquations[0] = x => -beta*x[1]*x[0];

                double [] initVal = {1,1.27e-6, 0}; //set inital values
                double time = 200; //set length of time
                double h = 1; //set timestep values
                RungeKutta rk = new RungeKutta(diffEquations,initVal,time,h); //pass all into Runge Kutta solver
                rk.Solve4thOrder(); //solve the system using 4th order
                String header = "Time,Suseptipal,Infected,Recovered"; //define file header
                rk.SaveResultsToCSV("results",header); //save results to a csv called results
            }catch(Exception e){
                Console.WriteLine(e.Message); //show error message to user
            }
        }
    }
    /*
    Implimintation of 4th order Runge Kutta Method
    */
    class RungeKutta{
        private functions[] equations = null; //set equations array to null so can be any size
        private functions[] Equations{ //getter setter for equations array
            get{return equations;}
            set{
                if(value.Length < 0){ //check array actually hase at least 1 eqn in it
                    throw new Exception("There must be at least 1 Differntial Equation");
                }
                equations = value;
            }
        }
        private double time = 0; //create time variable for length of time sim will run
        private double Time{
            get{return time;}
            set{
                if(value <= 0){ //check given time is positive
                    throw new Exception("Simulation time must be greater than 0");
                }
                time = value;
            }
        }
        private double h = 0; //set timestep
        private double H{
            get{return h;}
            set{
                if(value <= 0){ //check it is positive
                    throw new Exception("Timestep must be positive");
                }
                if(value > 5 ){ //if the value is greater than
                    Console.ForegroundColor = ConsoleColor.Yellow; //give warning that results may be low resulution
                    Console.WriteLine("WARNING: Large Timestep can lead to low resulution results");
                    Console.ForegroundColor = ConsoleColor.White;
                }
                h = value;
            }
        }
        private double Xj = 0; //set current timestep to 0
        private Vector[] solution = null; //define array of vectors to hold results
        public Vector[] Solution{get{return solution;}} //getter so user can use results
        private bool solved = false; //set solved to false
        private Vector K1,K2,K3,K4,initVals = null; //setup K vectors and inital value vector 
        private Vector InitVals{ // getter setter for inital values vector
            get{return initVals;}
            set{
                if(value.Size != Equations.Length){ //check there are the same number initali values as equations
                    throw new Exception("Must have same number of intial values as number of equations");
                }
                initVals = value;
            }
        }

        /*
        Constructor for RungeKutta solver
        */
        public RungeKutta(functions[] equations,double[] initVals,double time,double h){
            Equations = equations; //define equations to solve
            InitVals = new Vector(initVals); //define the inital values
            Time = time; //set time
            H = h; //set timestep
        }
        /*
        Solves the given system of equations using 4th Order Runge Kutta

        INPUT: | OUTPUT:
        Null   | NULL
        */
        public void Solve4thOrder(){
            Vector Yj = new Vector(equations.Length); //number of solutions will be number of equations
            Yj = initVals; //set inital values
            //define size of solution
            solution = new Vector[(int)(time/h)+1]; //plus 1 to garentee enough space casue dividion may not be whole number
            int itr = 0; //set iterations to zero
            
            while(Xj <= time){ //while having not reaching end time
                // First vector to save is the inital values
                Vector toSave = Yj.AddToFront(Xj); //add the timestep to the front of output
                solution[itr] = toSave; //add to solution array

                Vector K1 = K1Calc(Yj); //get K1 = f(Xj,Yj)
                Vector K2 = K2Calc(Yj); //get K2 = f(Xj+(h/2),Yj+h(K1/2))
                Vector K3 = K3Calc(Yj); //get K3 = f(Xj+(h/2),Yj+h(K2/2))
                Vector K4 = K4Calc(Yj); //get K4 = f(Xj+h,Yj+h*(K3))
                Yj = Yj + (1/6.0)*(h*(K1+(2.0*K2)+(2.0*K3)+K4)); //calulaute vector according to runge-kutta
                Xj = Xj+h; //inciment time
                itr++; //increase iteration
            }
            solved = true; //after loop the equations have been solved
        }

        /*
        Saves the solution of the system into a CSV file 

        INPUT:                                                           | OUTPUT:
        path - String of the location where the file should be created   | NULL
        */
        public void SaveResultsToCSV(string path,string header){
            StreamWriter s = new StreamWriter(path); //make new stream writer object
            if(solved == false){ //check if equation system has been solved yet
                throw new Exception("Cannot save results of an unsolved system");
            }
            if(!path.Contains(".csv")){ //if .csv not specified
                path = path+".csv"; //add the extension
            }

            string[] head = header.Split(','); //split header into string array
            if(head.Length != solution[1].Size){ //check that correct size header has been given
                Console.WriteLine("--------------------------------------------------------------------------");
                Console.WriteLine("Incorrect Number of Column Headers: File Will be saved without header line");
                Console.WriteLine("--------------------------------------------------------------------------");
            }else{
                s.WriteLine(header);
            }
            for(int i = 0; i < solution.Length; i++){ //for every vector in the solutions array
                s.WriteLine(solution[i]); //write to the file
            }
            Console.WriteLine("Saved File to: "+path); //let user know results have been saved
            s.Close(); //close file
        }

        /*
        Calculates the k1 vector for 4th order runge kutta
        K1 = f(Xj,Yj)
    
        INPUT:                                                           | OUTPUT:
        path - String of the location where the file should be created   | NULL
        */
        private Vector K1Calc(Vector Yj){ //works because there will always be the same number of equations as values
            K1 = new Vector(Yj.Size); //define Vector of size of the dependants pass thought to it 
            for(int i = 0; i < Yj.Size;i++){ //for each value
                K1[i] = equations[i](Xj,Yj); //each value in k corresponds to each eqn output
            }
            return K1; //return result
        } 

        /*
        Calculates the k2 vector for 4th order runge kutta 
        K2 = f(Xj+h/2,Yj+h(K1/2))

        INPUT:                                                           | OUTPUT:
        path - String of the location where the file should be created   | NULL
        */
        private Vector K2Calc(Vector Yj){
            K2= new Vector(Yj.Size); ////define Vector of size of the dependants pass thought to it 
            for(int i = 0; i < Yj.Size;i++){ //for each element in SIR
                K2[i] = equations[i](Xj+(h/2.0),Yj+(h*(K1/2.0))); //each value in k corresponds to each eqn output
            }
            return K2;
        } 

        /*
        Calculates the k3 vector for 4th order runge kutta 
        K3 = f(Xj+h/2,Yj+h(K2/2))

        INPUT:                                                           | OUTPUT:
        path - String of the location where the file should be created   | NULL
        */
        private Vector K3Calc(Vector Yj){
            K3 = new Vector(Yj.Size); //define size of Vector
            for(int i = 0; i < Yj.Size;i++){ //for each element in SIR
                K3[i] = equations[i](Xj+(h/2.0),Yj+(h*(K2/2.0))); //each value in k corresponds to each eqn output
            }
            return K3;
        } 

        /*
        Calculates the k4 vector for 4th order runge kutta 
        K4 = f(Xj+h,Yj+h*(K3))

        INPUT:                                                           | OUTPUT:
        path - String of the location where the file should be created   | NULL
        */
        public Vector K4Calc(Vector Yj){
            K4 = new Vector(Yj.Size); //define size of Vector
            for(int i = 0; i < Yj.Size;i++){ //for each element in SIR
                K4[i] = equations[i](Xj+h,Yj+(h*K3)); //each value in k corresponds to each eqn output
            }
            return K4;
        } 
    }
    class SIR_Model{
        private double d = 0; //defineing d value
        private double D{ //getter setter for D
            get{return D;}
            set{
                if(value <= 0 ){ //check to see if positive
                    throw new Exception("People must stay infected for positive length of time");
                }
                d = value;
            }
        }
        private double r0 = 0; //defining R0 value
        private double R0{
            get{return r0;}
            set{
                if(value <= 0 ){ //check to see if positive
                    throw new Exception("R0 must be positive");
                }
                r0 = value;
            }
        }
        private double beta = 0; //defining beta
        public double Beta{ //have public so user can overwrite value if they want to
            get{return beta;}
            set{
                if(value <= 0 ){ //check if positive
                    throw new Exception("beta must be positive");
                }
                beta = value;
            }
        }
        private double gamma = 0; //defining gamma
        public double Gamma{ //have public so user can overwrite value if they want to
            get{return gamma;}
            set{
                if(value <= 0 ){ //check if gamma is positive
                    throw new Exception("gamma must be positive");
                }
                gamma = value;
            }
        }  
        /*
        Constructor for SIR model
        */
        public SIR_Model(double D,double r0){
            this.D = D;
            gamma = 1.0/D;
            R0 = r0;
            beta = R0*gamma;
        }
        /*
        SIR model Equations
        values[0] = S
        values[1] = I
        values[2] = R

        INPUT:                           | OUTPUT:
        t = point in time                | single value from equation
        values = vector of SIR values
        While t is not striclty nessissary for theses equations it is good practise as other systems may not be autonomous
        */
        public double DSdt(double t,Vector values){
            return -beta*values[1]*values[0];///deps.sum();
        }

        public double DIdt(double t,Vector values){
            return ((beta*values[1]*values[0])-(gamma*values[1]));
        }

        public double DRdt(double t,Vector values){
            return gamma*values[1];
        }
    }

    /*
    Modified Vector class with redundant unused methods removed
    */
    public class Vector{
        private double [] data = null; //define data
        private int size;//define size
        public int Size{
            get{return size;}
            set{
                if(value < 1){
                    throw new Exception("Cannot make vector less than size 1");
                }
                size = value;
            }
        } 
        //getter setter for values
        public double this[int i]{
            get{return data[i];}
            set { data[i] = value; }
        }

        /*
        Initalise Two different ways of 
        */
        public Vector(int size){
            Size = size;
            data = new double[size];
            for(int i = 0; i < size;i++){
                data[i] = 0;
            }
        }
        public Vector(double[] data){
            Size = data.Length;
            this.data = new double [data.Length];
            for(int i =0; i < data.Length;i++){
                this.data[i] = data[i];
            }
        }

        /*
        Overloads * operator and Multiplies a vector by a double value element wise
        
        INPUT:                                                    | OUTPUT: 
        left = double that will be multiplied onto each element   | tmp = output vector with multiplied values
        right = Vector of doubles                                 | e.g left = a
                                                                        right = {b,c,d}
                                                                        left*right = {a*b,a*c,a*d}
        */ 
        public static Vector operator*(double left,Vector right){
            Vector tmp = new Vector(right.Size);
            for(int i = 0; i < right.Size; i++){
                tmp.data[i] = right.data[i]*left;
                //tmp.OverWrite(i,);
            }
            return tmp;
        }

        /*
        Overloads + operator and adds 2 vectors element wise
        
        INPUT:                     | OUTPUT: 
        left = Vector of doubles   | tmp = output vector with added values
        right = Vector of doubles  | e.g left = {a,b,c}
                                         right = {d,e,f}
                                         left+right = {a+d,b+e,c+f}
        */ 
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
        Overloads / operator and divids each element of a vector by a double
        
        INPUT:                          | OUTPUT: 
        left = Vector of doubles        | tmp = output vector with added values
        right = doubles to divide with  | e.g left = {a,b,c}
                                              right = d
                                              left/right = {a/d,b/d,c/d}
        */  
        public static Vector operator/(Vector left, double right){
            Vector tmp = new Vector(left.Size);
            for(int i = 0;i < tmp.Size;i++){
                tmp.data[i] = left.data[i] / right;
            }
            return tmp;
        } 

        /*
        Adds a double onto the front of a vector
        
        INPUT:                              | OUTPUT: 
        number = double to add to front     | tmp = Vector with value added
                                            | e.g Vector = {a,b,c}
                                                  Vector.AddToFront(d) = {d,a,b,c}
        */ 
        public Vector AddToFront(double number){
            Vector tmp = new Vector(Size+1); //make new vector one size bigger than current
            tmp[0] = number; //put number in first element
            for(int i = 1; i < tmp.Size;i++){ //loop through rest for current values
                tmp[i] = data[i-1]; //add to tmp
            }
            return tmp; //return tmp
        }

        /*
        Generates string version of vector overloading ToString() method
        
        INPUT:  | OUTPUT: 
        NULL    | tmp = string representaion on Vector
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

    }
}
