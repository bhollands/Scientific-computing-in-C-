/*
Author: Bernard Hollands

Gauss-Jacobi Solver in C# written for Assignment 4 of AM6007

Program.cs File contains 4 Classes

Progam - Driver code with main method identical to given example in breif
Matrix - Impliments a sqaure matrix where data is held in 2D Array
Vector - Impliments vectors where data is held in 1D Array
Linsolve -  Impliments Gauss-Jacobi method for solving linear equations
*/


using System;
namespace Gauss_Jacobi
{
    class Program
    {
        static void Main(string[] args)
        {
            /*
            Exact code from breif
            */
            try
            {
                Matrix m = new Matrix(4); //make 4x4 matrix
                Vector b = new Vector(4); //make 4x1 vector
                b[0] = 54.5; b[1] = -14; b[2] = 12.5; b[3] = -21; //define vector elements
                m[0, 0] = 9; m[0, 1] = -2; m[0, 2] = 3; m[0, 3] = 2; //define Matrix elements
                m[1, 0] = 2; m[1, 1] = 8; m[1, 2] = -2; m[1, 3] = 3;
                m[2, 0] = -3; m[2, 1] = 2; m[2, 2] = 11; m[2, 3] = -4;
                m[3, 0] = -2; m[3, 1] = 3; m[3, 2] = 2; m[3, 3] = 10;
                Console.WriteLine("The Matrix m is{0}", m); //output matrix to console
                Console.WriteLine("The Vector b is: {0}", b); //output vector to console
                Linsolve l = new Linsolve(100); // create a linear solver with 100 iterations
                Vector ans = l.Solve(m, b); //solve m and b
                Console.WriteLine("The Solution to mx = b is {0}", ans); //output solution to console
            }
            catch (Exception e) //catch any excetions
            {
                Console.WriteLine("Error encountered: {0}", e.Message); //display error if one occured
            }
            Console.ReadLine();

        }
    }

    //Linear solver class 
    class Linsolve
    {
        private int maxItr; //max iterations allowed
        // getter setter for max iteration
        public int MaxItr 
                        {
                            get{return maxItr;}
                            set
                            {
                                if(value < 0) //if max iterations is negative
                                {
                                    throw new Exception("Cannot set negative max iterations");
                                }
                                maxItr = value;
                            }
                        }
        //Construtor to set max iterations
        public Linsolve(int maxItr)
        {
            MaxItr = maxItr;
        }
        // default constructor sets max iterations to 100
        public Linsolve()
        {
            MaxItr = 100;
        }
        //define matrices that need to be calculated
        private Matrix T = null;
        private Vector c = null;
        private Matrix D = null;
        private Matrix U = null;
        private Matrix L = null;

        /*
        Solves a system of linear equations using Gauss-Jacobi Method

        INPUT:                                   | OUTPUT:  
        A = Input matrix from equaions           | xk1 = Vector of the found unknowns
        b = Input Vector from equation solutions | 
        */
        public Vector Solve(Matrix A, Vector b)
        {
            FindDeps(A, b); //Generate the required matices and Vectors for Gauss-Jacobi
            int itr = 0; //iteration counter
            Vector xk = new Vector(b.Size); // Vector to hold current value
            Vector xk1 = new Vector(b.Size); //Vector to hold next iteration
            for(int i = 0; i <= MaxItr; i++)
            {
                xk1 = T * xk + c; //calculate next step
                double error = (xk1 - xk).norm() / xk.norm(); //calclate checking number
                if (error < 10e-7) //if error is small enough
                {
                    return xk1.Round(3);//return the solution rounded to 3 decimal places
                }
                xk = xk1; //set ccurrent next step to current value
                itr++; //itercreate iteration counter
            }
            // if loop completed without convergence throw exception
            throw new Exception("Unable to Converge after " + MaxItr + " Iterations"); //Let user know
        }
        /*
        Finds the required matrices and vectors

        INPUT:                                   | OUTPUT:  
        A = Input matrix from equaions           | void
        b = Input Vector from equation solutions | 
        */
        private void FindDeps(Matrix A, Vector b)
        {
            D = !(~A); //inverse the diagonal
            L = -Matrix.BelowDi(A); //minus below diagonal elements
            U = -Matrix.AboveDi(A); //minus the elemets above the diagonal
            T = D * (L + U); //calculate T as per equation 
            c = D * b; //calculate c as Per equation
        }
    }
    /*
    Matrix class using 2D array
    */
    public class Matrix
    {

        //make better so can get entire rows and columns

        private int size; //size, since nxn size is both length and width
        private double[,] matrix = null; //define null 2D matrix to hold data
        // getter setter for size
        public int Size
        {
            get { return size; }
            set
            {
                if (value < 1) //if trying to set size less than 1
                {
                    throw new Exception("Cannot make Matrix of less than size 1");
                }
                size = value; //set size
                matrix = new double[size, size];
            }
        } //return size as a vector

        //indexer for individual values for the matrix
        public double this[int row, int col]
        {
            get
            {
                CheckRange(row,col);
                return matrix[row, col];
            }
            set
            {
                CheckRange(row,col);
                matrix[row, col] = value;
            }
        }
        private void CheckRange(int row, int col)
        {
            if (row < 0 || col < 0) //same checks above
                {
                    throw new Exception("Cannot set negative matrix index");
                }
                else if (row >= Size || col >= Size)
                {
                    throw new Exception("Cannot set matrix index greater than Matrix size");
                }
        }
        // overloading of the indexing function to set entire rows and columns from vectors
        private string validRequests = "r,c"; //string of value chars
        //r - reuqets a row
        //c - requets a column 
        public Vector this[int index, char a]
        {
            get
            {
                Vector tmp = new Vector(Size); //make new vector to return
                if(!validRequests.Contains(a)){throw new Exception(a + " is not a a valid matrix request");} //check it is a valid request
                if(a == 'r'){tmp = GetRow(index);} //get the row specified
                if(a == 'c'){tmp = GetColumn(index);} //get the column specified
                return tmp; //return the vector
            }
            set
            {
                if(!validRequests.Contains(a)){throw new Exception(a + " is not a a valid matrix request"); }//check it is a valid request
                if(a == 'r'){SetRow(index, value);} //set row specified
                if(a == 'c'){SetColumn(index, value);} //set column specified
            }
        }

        //Constructor to set matrix size
        public Matrix(int size)
        {
            Size = size; //set size
            InitMatrix(); //initalise the matrix with all zeros
        }
        // defult constructor to set size 3 matrix
        public Matrix()
        {
            Size = 3;
            InitMatrix(); //initalise the matrix with all zeros

        }
        private void InitMatrix(){
            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < Size; j++)
                {
                    matrix[i, j] = 0;
                }
            }
        }

        /*
        Gets an entire column form matrix and return it as a vector

        INPUT:                                  | OUTPUT:  
        index = index of column to make vector   | tmp = vector of value from the column at index
        */
        private Vector GetColumn(int index)
        {
            if(index < 0 || index > Size) //check index is in valid range
            {  
                throw new Exception("Cannot get column at index "+ index +"in matrix of size"+Size);
            }
            Vector tmp = new Vector(Size); //make vector to return
            for (int i = 0; i < Size; i++)//for the size of the matrix
            {
                tmp[i] = matrix[i, index]; //get the value from the column at each value 1
            }
            return tmp;
        }
        /*
        Gets an entire row form matrix and return it as a vector

        INPUT:                                  | OUTPUT:  
        index = index of row to make vector   | tmp = vector of values from the row at index
        */
        private Vector GetRow(int index)
        {
            if(index < 0 || index > Size) //check index is valid
            {  
                throw new Exception("Cannot get row at index "+ index +"in matrix of size"+Size);
            }
            Vector tmp = new Vector(Size);//make vector to return
            for (int i = 0; i < Size; i++)
            {
                tmp[i] = matrix[index, i]; //get the value from the row at each value 1
            }
            return tmp; //return the vector
        }
        /*
        Sets an entire row in a matrix from a vector

        INPUT:                      | OUTPUT:  
        index = index of row to set | void
        */
        private void SetRow(int index, Vector value)
        {
            CheckIndex(index,value);
            for (int i = 0; i < Size; i++) //for the size
            {
                matrix[index, i] = value[i]; //set each value in the index row
            }
        }

        /*
        Sets an entire column in a matrix from a vector

        INPUT:                                  | OUTPUT:  
        index = index of column to ake vector   | void
        */
        private void SetColumn(int index, Vector value)
        {
            CheckIndex(index,value);
            for (int i = 0; i < Size; i++)
            {
                matrix[i, index] = value[i]; //set each value in the index columns
            }
        }
        private void CheckIndex(int index,Vector Value)
        {
            if(index < 0 || index > Size) //check for valid index
            {  
                throw new Exception("Cannot set column at index "+ index +"in matrix of size"+Size);
            }
        }
        //Override of ToString() method
        public override string ToString()
        {
            string tmp = "\n["; //set string to hold matrix representation
            
            for (int i = 0; i < Size; i++) //for the size of the matrix (rows)
            {
                string delimiter = ",";// set delimiter between numbers
                for (int j = 0; j < Size; j++) //for the size of the matrix (rows)
                {
                    if (j == Size - 1)//if at end of columns
                    { 
                        delimiter = "]\n["; //put end cap and drop a line with opening cap [
                        if (i == Size - 1) //if at final end
                        {
                            delimiter = "]\n"; //put end cap a drop line
                        }
                    }
                    tmp = tmp + matrix[i, j] + delimiter;// add to string prespresentation
                }
            }
            return tmp;//return the strings
        }

        //--- Static Functions ------------------

        /*
        Returns the top left to bottom right diagonal elements of a matrix

        INPUT:                                  | OUTPUT:  
        A = Matrix to find diagonal elements    | tmpMat = Matrix with only diagonal elements preserved
                                                | e.g A = [a,b,c] tmpMat = [a,0,0] 
                                                          [d,f,g]          [0,f,0]
                                                          [h,i,j]          [0,0,j]
        */
        public static Matrix operator ~(Matrix A)
        {
            Vector tmpVec = null;
            Matrix tmpMat = new Matrix(A.Size);
            for (int i = 0; i < A.Size; i++) //for the size of A
            {
                tmpVec = A[i,'r']; //get the vector
                tmpMat[i,'r'] = tmpVec.Preserve(i); //make all zero but i
            }
            return tmpMat;
        }
        /*
        Returns the all elements above the diagonal line

        INPUT:                                        | OUTPUT:  
        A = Matrix to find above diagonal elements    | tmpMat = Matrix with only above diagonal elements preserved
                                                      | e.g A = [a,b,c] tmpMat = [0,b,c] 
                                                                [d,f,g]          [0,0,g]
                                                                [h,i,j]          [0,0,0]
        */
        public static Matrix AboveDi(Matrix A)
        {
            Matrix tmpMat = new Matrix(A.Size);
            for (int i = 0; i < A.Size; i++) //for size of A
            {
                Vector tmpVec = A[i,'r']; //get the current vector
                tmpMat[i,'r'] = tmpVec.ZerosRange(0, i); //make elements from 0 to i Zero
            }
            return tmpMat; //return mat
        }
        /*
        Returns the all elements below the diagonal line

        INPUT:                                        | OUTPUT:  
        A = Matrix to find below diagonal elements    | tmpMat = Matrix with only below diagonal elements preserved
                                                      | e.g A = [a,b,c] tmpMat = [0,0,0] 
                                                                [d,f,g]          [d,0,0]
                                                                [h,i,j]          [h,i,0]
        */
        public static Matrix BelowDi(Matrix A)
        {
            Matrix tmpMat = new Matrix(A.Size);
            for (int i = 0; i < A.Size; i++)
            {
                Vector tmpVec = A[i,'r']; //get the vector
                tmpMat[i,'r'] = tmpVec.ZerosRange(i, A.Size); //make everything in range zero
            }
            return tmpMat;
        }

        /*
        Returns Matrix where all elements have been inversed

        INPUT:                         | OUTPUT:  
        A = Matrix to find inveres of  | tmpMat = Matrix with all elements inversed
                                       | e.g A = [a,b,c] tmpMat = [1/a,1/b,1/c] 
                                                 [d,f,g]          [1/d,1/f,1/g]
                                                 [h,i,j]          [1/h,1/i,1/j]
        */
        public static Matrix operator !(Matrix A)
        {
            Vector tmp = null;
            Matrix tmpMat = new Matrix(A.Size);
            for (int i = 0; i < A.Size; i++)
            {
                tmp = A[i,'r']; //get the vector
                tmp = tmp ^ ((1 / tmp) / tmp); //inverse | ^ is element wise multiiplication
                tmpMat[i,'r'] = tmp;
            }
            return tmpMat;
        }
        /*
        Sums 2 Matricses element wise

        INPUT:                 | OUTPUT:  
        lhs = Matrix to add    | tmp = Matrix with the element wise sum
        rhs = Matrix to add    | e.g lhs = [a,b,c] rhs = [1,2,3] tmp = [a+1,b+2,c+3]
                                           [d,f,g]       [4,5,6]       [d+4,f+5,g+6]
                                           [h,i,j]       [7,8,9]       [h+7,8+i,9+j]
        */
        public static Matrix operator +(Matrix lhs, Matrix rhs)
        {
            if (lhs.Size != rhs.Size)
            {
                throw new Exception("Unable to Sum Matrices of size " + lhs.Size + " and " + rhs.Size);
            }
            Matrix tmp = new Matrix(rhs.Size);
            for (int i = 0; i < tmp.Size; i++)
            {
                tmp[i,'r'] = lhs[i,'r'] + rhs[i,'r'];
            }
            return tmp;
        }
        /*
        Makes all elements of a matrrix negtive

        INPUT:                         | OUTPUT:  
        lhs = Matrix to make negative  | tmpMat = Matrix with all elements inversed
                                       | e.g A = [a,b,c] tmpMat = [-a,-b,-c] 
                                                 [d,f,g]          [-d,-f,-g]
                                                 [h,i,j]          [-h,-i,-j]
        */
        public static Matrix operator -(Matrix lhs)
        {
            Matrix tmp = new Matrix(lhs.Size);
            for (int i = 0; i < tmp.Size; i++)
            {
                tmp[i,'r'] = -lhs[i,'r'];
            }
            return tmp;
        }

        /*
        Cross multiplies 2 matices

        INPUT:                      | OUTPUT:  
        lhs = Matrix to multiply    | tmp = Matrix with the element wise sum
        rhs = Matrix to multiply    | e.g lhs = [a,b] rhs = [1,2] tmp = [(a*1)+(b*4),(a*2)+(b*5)]
                                                [d,f]       [4,5]       [(d*1)+(f*4),(d*2)+(f*5)]                                  
        */
        public static Matrix operator *(Matrix lhs, Matrix rhs)
        {
            if (lhs.Size != rhs.Size) //check operation is valid
            {
                throw new Exception("Unable to Multiply Matrices of size " + lhs.Size + " and " + rhs.Size);
            }
            Matrix tmp = new Matrix(rhs.Size); // make tmp matrix
            double dot = 0;// dot product result
            for (int i = 0; i < tmp.Size; i++) //for the size of matrix (rows)
            {
                for (int j = 0; j < tmp.Size; j++) //for the size of matrix (columns)
                {
                    dot = lhs[i,'r'] * rhs[j,'c'];//for row i in lhs find its dot prodcut with with columns of rhs
                    tmp[i,j] = dot; //store in output matrix
                }
            }
            return tmp;
        }
        /*
        Cross multiplies matrix and vector

        INPUT:                      | OUTPUT:  
        lhs = Matrix to multiply    | tmp = Matrix with the element wise sum
        rhs = Matrix to multiply    | e.g lhs = [a,b] rhs = [1,2] tmp = [(1*a)+(2*b),(1*c)+(2*d)]
                                                [d,f]                                                  
        */
        public static Vector operator *(Matrix lhs, Vector rhs)
        {
            if (lhs.Size != rhs.Size) //check operation is valid
            {
                throw new Exception("Unable to Multiply Matrices of size " + lhs.Size + " and Vector of Size" + rhs.Size);
            }
            Vector tmp = new Vector(rhs.Size); //make output vector
            for (int i = 0; i < tmp.Size; i++) //for the size of output vector
            {
                tmp[i] = lhs[i,'r'] * rhs; //dot product each matrux row with vector
            }
            return tmp;
        }

    }

    //Vector class
    public class Vector
    {
        private double[] data = null; //array to hold data
        //getter setter for size
        public int Size
        {
            get { return data.Length; }
            set
            {
                if (value < 1)
                {
                    throw new Exception("Cannot make Vector of less than size 1");
                }
                data = new double[value];
            }
        } //return size as a vector
        
        //indexing function
        public double this[int i]
        {
            get
            {
                CheckRange(i);
                return data[i];
            }
            set
            {
                CheckRange(i);
                data[i] = value;
            }
        }

        private void CheckRange(int i)
        {
            if (i < 0 || i > Size)
                {
                    throw new Exception("Cannot accsess vector positions outside range 0->" + Size);
                }
        }
        //default constructor
        public Vector()
        {
            Size = 3;//set size
            for (int i = 0; i < 3; i++)
            {
                data[i] = 0;
            }
        }
        // initalise vector of certin size
        public Vector(int size)
        {
            Size = size; //set size
            for (int i = 0; i < Size; i++)
            {
                data[i] = 0; //initalise with zeros
            }
        }
        //initalise vector from array
        public Vector(double[] data)
        {
            Size = data.Length; //set size
            for (int i = 0; i < data.Length; i++) //for each element
            {
                this.data[i] = data[i]; //copy the information over
            }
        }


        /*
        Find Scalar/Dot product of 2 vectors

        INPUT:                      | OUTPUT:  
        lhs = Vector to multiply    | tmp = double with the element wise sum
        rhs = Vector to multiply    | e.g lhs = [a,b] rhs = [1,2] tmp = (1*a)+(2*b)                                                 
        */
        public static double operator *(Vector lhs, Vector rhs)
        {
            CheckSize(rhs,lhs);
            double tmp = 0;
            for (int i = 0; i < lhs.Size; i++) //for the size of the vectors
            {
                tmp += lhs.data[i] * rhs.data[i]; //sum the product of each element in index i
            }
            return tmp;
        }

        /*
        Calulates elemens wise multiplication of 2 vectors

        INPUT:                      | OUTPUT:  
        lhs = Vector to multiply    | tmp = Vector with the element wise multiplications
        rhs = Vector to multiply    | e.g lhs = [a,b] rhs = [1,2] tmp = [1*a,2*b]                                                 
        */
        public static Vector operator ^(Vector left, Vector right)
        {
            CheckSize(left,right);
            Vector tmp = new Vector(left.Size); ;
            for (int i = 0; i < left.Size; i++) //for each element
            {
                tmp[i] = left.data[i] * right.data[i]; //multiply each element and put into output vector
            }
            return tmp;
        }
        /*
        Outputs matrix where each element is the negative version

        INPUT:                      | OUTPUT:  
        left = Vector to multiply   | tmp = Vector with the negative elements
                                    | e.g left = [a,b,c] tmp = [-a,-b,-c]
        */
        public static Vector operator -(Vector left)
        {
            Vector tmp = new Vector(left.Size);
            for (int i = 0; i < left.Size; i++)
            {
                tmp.data[i] = -left.data[i];//make minus
                //tmp.OverWrite(i,);
            }
            return tmp;
        }

        /*
        Outputs matrix where each element is the element wise sub-traction between vectors

        INPUT:                           | OUTPUT:  
        a = Vector on lhs to subtract    | tmp = Vector with the sub tracted elements
        b = Vector on rhs in subtraction | e.g a = [a,b,c] b = [d,e,f], tmp = [a-d,b-e,c-f]
        */
        public static Vector operator -(Vector a, Vector b)
        {
            CheckSize(a,b);
            Vector tmp = new Vector(a.Size);
            for (int i = 0; i < a.Size; i++) //for every element
            {
                tmp.data[i] = a[i] - b[i]; //subtract each element from each other
            }
            return tmp;
        }

        /*
        Returns the abolute value of the Vector 

        INPUT:  | OUTPUT:  
        void    | result = absolute value
                | e.g [a,b,c].norm() = sqrt(a^2 + b^2 + c^2)
        */
        public double norm()
        {
            double result = 0;
            for (int i = 0; i < Size; i++)//for each element
            {
                result = result + (data[i] * data[i]); //sum the squres of each element 
            }
            result = Math.Sqrt(result); //find the squre root of the sum
            return result;
        }

        /*
        Returns the element wise sum of 2 vectors 

        INPUT:               | OUTPUT:  
        a = Vector to add    | tmp = Vector with element wise sums
        b = Vector to add    | e.g e.g a = [a,b,c] b = [d,e,f], tmp = [a+d,b+e,c+f]
        */
        public static Vector operator +(Vector a, Vector b)
        {
            CheckSize(a,b);
            Vector tmp = new Vector(a.Size);
            for (int i = 0; i < tmp.Size; i++) //for every element
            {
                tmp.data[i] = a.data[i] + b.data[i]; //sum each element
            }
            return tmp;
        }

        /*
        INPUT:                          | OUTPUT:  
        left = double on the numerator  | tmp = Vector with division results
        right = Vector to be divided    | e.g e.g left = 1 right = [d,e,f], tmp = [1/d,1/e,1/f]
        */
        public static Vector operator /(double left, Vector right)
        {
            if(left == 0) //check for zero divide
            {   
                throw new DivideByZeroException("Cannot divide by zero");
            }
            Vector tmp = new Vector(right.Size);
            for (int i = 0; i < tmp.Size; i++) //for the size of the vector
            {
                tmp.data[i] = left / right.data[i]; //divide left by every element
            }
            return tmp;
        }
        /*
        Divdes 2 Vector element wise

        INPUT:                          | OUTPUT:  
        left = Vector on the numerator  | tmp = Vector with division results
        right = Vector on denominator   | e.g e.g left = [a,b,c] right = [d,e,f], tmp = [a/d,b/e,c/f]
        */
        public static Vector operator /(Vector left, Vector right)
        {
            CheckSize(left,right); //check size
            Vector tmp = new Vector(right.Size); //make new vector
            for (int i = 0; i < tmp.Size; i++) //for size of vector
            {
                if (right[i] != 0) //if not equal to zero
                {
                    tmp.data[i] = left.data[i] / right.data[i]; //do the dividion
                }
            }
            return tmp;
        }

        /*
        Checks 2 Vectors are the same size if now throws exception

        INPUT:                  | OUTPUT:  
        left = Vector to check  | void
        right = Vector to check |
        */
        private static void CheckSize(Vector left, Vector right)
        {
            if (left.Size != right.Size) // check same size
            {
                throw new Exception("Cannot add Vectors of size " + right.Size + " and " + right.Size); // if not throw exception
            }
        }
        /*
        Makes all elements in a vector zero but one

        INPUT:                | OUTPUT:  
        index = index to save | tmp = Vector with preserved element
        */
        public Vector Preserve(int index)
        {
            CheckRange(index); //check index is in correct range
            for (int i = 0; i < Size; i++)
            {
                if (i != index) // if not the specifed index
                {
                    data[i] = 0; //set to zero
                }
            }
            Vector tmp = new Vector(data);
            return tmp;
        }

        /*
        Makes all emements in the range zero

        INPUT:              | OUTPUT:  
        start = start index | tmp = Vector with Range zeroed
        end = end index
        */
        public Vector ZerosRange(int start, int end)
        {
            if(start > end || start < 0 || end > Size)
            { 
                throw new Exception("Invalid range to make zero");
            }
            //make new vector
            for (int i = 0; i < Size; i++) //for the size of the vector
            {
                if (i >= start && i <= end) //if in range
                {
                    data[i] = 0; //set to zero
                }
            }
            Vector tmp = new Vector(data);
            return tmp;
        }
        /*
        Rounds all elements in Vector to a specifed precision

        INPUT:                                            | OUTPUT:  
        precision = Number of decimal places to round to  | tmp = Vector with rounded values
        end = end index
        */
        public Vector Round(int precision)
        {
            if(precision < 0){ //check valid range
                throw new Exception("Cannot have negative number of decimal places");
            }
            for (int i = 0; i < data.Length; i++) //for the size of the vector
            {
                data[i] = Math.Round(data[i], precision); //round element to precision
            }
            Vector tmp = new Vector(data);
            return tmp;
        }
        /*
        Overloads to string class to print vector
        */
        public override string ToString()
        {
            string tmp = "["; //set temp string
            for (int i = 0; i < Size; i++)
            { //for the size of the vector
                tmp += string.Format("{0:0.000}", data[i]); //for each element turn into string and add to tmp
                if (i != Size - 1)
                { //if not at end position
                    tmp += ","; //add a comman
                }  //This is avoid a comma after the final number on printing
            }
            return tmp + "]"; //return the string
        }
    }
}
