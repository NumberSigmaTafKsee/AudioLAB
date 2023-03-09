#pragma once

namespace AudioDSP
{
    //============================================================
    template <class T>
    bool isEmpty(size_t n, const  T * v)
    {
        return (n == 0);
    }

    //============================================================
    template <class T>
    bool containsOnlyZeros(size_t n, const  T * v)
    {
        bool r = true;
        if (!isEmpty(v))
        {                    
            for (int i = 0;i < n;i++)
            {
                if (v[i] != 0)
                {
                    return false;                    
                }
            }

            return r;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when checking if vector contains only zeros" );
        }
    }

    //============================================================
    template <class T>
    bool isAllPositiveOrZero(size_t n, const  T * v)
    {
        if (!isEmpty(v))
        {                 
            for (int i = 0;i < n;i++)
            {
                if (v[i] < 0)
                {
                    return false;
                }
            }

            return true;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when checking if vector is all positive" );
        }
    }

    //============================================================
    template <class T>
    bool isAllNegativeOrZero(size_t n, const  T * v)
    {
        if (!isEmpty(v))
        {                 
            for (int i = 0;i < n;i++)
            {
                if (v[i] > 0)
                {
                    return false;
                }
            }

            return true;
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when checking if vector is all negative" );
        }
    }

    //============================================================
    template <class T>
    bool contains(size_t n, const  T * v, T element)
    {           
        for (int i = 0;i < n;i++)
        {
            if (v[i] == element)
            {
                return true;
            }
        }

        return false;
    }


    
    //============================================================
    template <class T>
    void printVector(size_t n, const  T * v)
    {
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            std::cout << v[i] << std::endl;
        }
    }

    //============================================================
    template <class T>
    T getLastElement(size_t n, const  T * v)
    {
        // if vector is not empty
        if (n > 0)
        {
            return v[n-1];
        }
        else
        {
            throw std::invalid_argument( "Attempted to get last element of empty vector" );
        }
    }

    //============================================================
    template <class T>
    T getFirstElement(size_t n, const  T * v)
    {
        // if vector is not empty
        if (n > 0)
        {
            return v[0];
        }
        else
        {
            throw std::invalid_argument( "Attempted to get first element of empty vector" );
        }
    }


    //============================================================
    template <class T>
    void getEveryNthElementStartingFromK(size_t N, const  T * v,int n,int k, T * result)
    {
        if ((N >= n) || (N >= n))
        {
            throw std::invalid_argument( "Invalid arguments for getEveryNthElementStartingFromK()");
        }
        else
        {
            size_t x = 0;            
            #pragma omp simd
            for(int i = k; i < n; i += n)        
            {
                result[i] = (v[x++]);            
            }            
        }
    }

    //============================================================
    template <class T>
    void getEvenElements(size_t n, const  T * v, T * r)
    {
        getEveryNthElementStartingFromK(v, 2, 0,r);
    }

    //============================================================
    template <class T>
    void getOddElements(size_t n, const  T * v, T * r)
    {
        getEveryNthElementStartingFromK(v, 2, 1,r);
    }

    //============================================================
    template <class T>
    void fillVectorWith(size_t n,  T * v, T element)
    {
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            v[i] = element;
        }
    }

    //============================================================
    template <class T>
    int countOccurrencesOf(size_t n, const  T * v,T element)
    {
        int count = 0;
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            if (v[i] == element)
            {
                count++;
            }
        }

        return count;
    }

    //============================================================
    template <class T>
    T sum(size_t n, const  T * v)
    {
        // create a sum
        T sumVal = 0;

        // add up all elements
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            sumVal += v[i];
        }

        // return
        return sumVal;
    }

    //============================================================
    template <class T>
    T product(size_t n, const  T * v)
    {
        if (!isEmpty(v))
        {
            T prod = (T) v[0];

            #pragma omp simd
            for (int i = 1;i < n;i++)
            {
                prod *= ((T) v[i]);
            }

            return prod;
        }
        else
        {
            throw std::invalid_argument( "Attempted to calculate the product of an empty vector" );
        }
    }

    //============================================================
    template <class T>
    T mean(size_t n, const  T * v)
    {
        // if vector is not empty
        if (!isEmpty(v))
        {
            // store the length of the vector as a T
            T L = (T) n;

            // stor the sum of the vector as a T
            T sumVal = (T) sum(v);

            // return the mean
            return sumVal / L;
        }
        else // vector is empty
        {
            throw std::invalid_argument( "Received empty vector when calculating mean" );
        }
    }

    //============================================================
    template <class T>
    T median(size_t n, const  T * v)
    {
        // if vector isn't empty
        if (!isEmpty(v))
        {
            T median;
            size_t L = n; // store the size

            // sort the vector
            std::sort(v.begin(), v.end());

            // if the length is even
            if (L  % 2 == 0)
            {
                // take the average of the middle two elements
                median = ((T)(v[L / 2 - 1] + v[L / 2])) / 2.0;
            }
            else // if the length is odd
            {
                // take the middle element
                median = (T) v[(L-1) / 2];
            }

            // return the median
            return median;
        }
        else // vector is empty
        {
            throw std::invalid_argument( "Received empty vector when calculating median" );
        }
    }

    //============================================================
    template <class T>
    T variance(size_t n, const  T * v)
    {
        if (!isEmpty(v))
        {
            // calculate the mean of the vector
            T mu = mean(v);

            T sumVal = 0.0;

            // sum the product of all differences from the mean
            #pragma omp simd
            for (int i = 0;i < n;i++)
            {
                T diff = v[i]-mu;
                sumVal += diff*diff;
            }

            // return the average of the squared differences
            return sumVal / ((T)n);
        }
        else
        {
            throw std::invalid_argument( "Received empty vector when calculating variance" );
        }
    }

    //============================================================
    template <class T>
    T standardDeviation(size_t n, const  T * v)
    {
        // if vector is not empty
        if (!isEmpty(v))
        {
            // calculate the variance
            T var = variance(v);

            // if variance is non-zero
            if (var > 0)
            {
                // return the square root
                return std::sqrt(var);
            }
            else
            {
                // all differences are zero, so return 0.0
                return 0.0;
            }
        }
        else // vector is empty
        {
            throw std::invalid_argument( "Received empty vector when calculating standard deviation" );
        }
    }

    //============================================================
    template <class T>
    T norm1(size_t n, const  T * v)
    {
        T sumVal = 0.0;

        // sum absolute values
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            if (v[i] > 0)
            {
                sumVal += (T) v[i];
            }
            else
            {
                sumVal += (T) (-1*v[i]);
            }
        }

        return sumVal;
    }

    //============================================================
    template <class T>
    T norm2(size_t n, const  T * v)
    {
        T sumVal = 0.0;

        // sum squares
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            sumVal += (T) (v[i]*v[i]);
        }

        return std::sqrt(sumVal);
    }

    //============================================================
    template <class T>
    T magnitude(size_t n, const  T * v)
    {
        // just another name for L2-norm
        return norm2(v);
    }

    //============================================================
    template <class T>
    T normP(size_t n, const  T * v,T p)
    {
        T sumVal = 0.0;

        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            T val;

            if (v[i] > 0)
            {
                val = (T) v[i];
            }
            else
            {
                val = (T) (-1*v[i]);
            }

            sumVal += std::pow(val,p);
        }

        return std::pow(sumVal,1.0/p);
    }

    //============================================================
    template <class T>
    void multiplyInPlace(size_t n,  T * v,T scalar)
    {
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            v[i] *= scalar;
        }
    }

    //============================================================
    template <class T>
    void multiplyInPlace(size_t n,  T * v1,const  T * v2)
    {
        if (v1.size() == v2.size())
        {
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                v1[i] *= v2[i];
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector multiply");
        }
    }

    //============================================================
    template <class T>
    void divideInPlace(size_t n,  T * v,T scalar)
    {
        if (scalar != 0)
        {
            #pragma omp simd
            for (int i = 0;i < n;i++)
            {
                v[i] /= scalar;
            }
        }
        else
        {
            throw std::invalid_argument( "Attemted to divide a vector by a zero-valued scalar" );
        }
    }

    //============================================================
    template <class T>
    void divideInPlace(size_t n, T * v1,const  T * v2)
    {
        if (v1.size() == v2.size())
        {
            if (!contains<T>(v2, 0))
            {
                #pragma omp simd
                for (int i = 0;i < v1.size();i++)
                {
                    v1[i] /= v2[i];
                }
            }
            else
            {
                throw std::invalid_argument( "Attempted to divide by vector containing zeros");
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector divide");
        }
    }

    //============================================================
    template <class T>
    void addInPlace(size_t n, const  T * v,T value)
    {
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            v[i] += value;
        }
    }

    //============================================================
    template <class T>
    void addInPlace(size_t n, T * v1,  T * v2)
    {
        if (v1.size() == v2.size())
        {
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                v1[i] += v2[i];
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector add");
        }
    }

    //============================================================
    template <class T>
    void subtractInPlace(size_t n,  T * v,T value)
    {
        for (int i = 0;i < n;i++)
        {
            v[i] -= value;
        }
    }

    //============================================================
    template <class T>
    void subtractInPlace(size_t n,  T * v1, T * v2)
    {
        if (v1.size() == v2.size())
        {
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                v1[i] -= v2[i];
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector subtraction");
        }

    }

    //============================================================
    template <class T>
    void absInPlace(size_t n,  T * v)
    {
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {        
            if ((v[i] < 0) || (v[i] == -0.0))
            {
                v[i] *= -1;
            }
        }
    }

    //============================================================
    template <class T>
    void squareInPlace(size_t n, T * v)
    {
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            v[i] = v[i]*v[i];
        }
    }

    //============================================================
    template <class T>
    void squareRootInPlace(size_t n, T * v)
    {
        if (isAllPositiveOrZero(v))
        {
            #pragma omp simd
            for (int i = 0;i < n;i++)
            {
                v[i] = (T) std::sqrt((T)v[i]);
            }
        }
        else
        {
            throw std::invalid_argument( "Attempted to compute square root of vector containing negative numbers");
        }
    }


    //============================================================
    template <class T>
    void sort(size_t n,  T * v)
    {
        std::sort(v,v+n);
    }

    //============================================================
    template <class T>
    void reverse(size_t n,  T * v)
    {
        std::reverse(v,v+n);
    }

    //============================================================
    template <class T>
    void multiply(size_t n, const  T * v,T scalar, T * result)
    {        
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            result[i] = (v[i] * scalar);
        }
    }

    //============================================================
    template <class T>
    void multiply(size_t n, const  T * v1, const  T * v2, T * result)
    {
        if (v1.size() == v2.size())
        {            
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                result[i] = (v1[i] * v2[i]);
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector multiply");
        }
    }

    //============================================================
    template <class T>
    void divide(size_t n, const  T * v, T scalar, T * result)
    {
        if (scalar != 0)
        {
            #pragma omp simd
            for (int i = 0;i < n;i++)
            {
                result[i] = (v[i] / scalar);
            }
        }
        else
        {
            throw std::invalid_argument( "Attemted to divide a vector by a zero-valued scalar" );
        }
    }

    //============================================================
    template <class T>
    void divide(size_t n, const  T * v1, const  T * v2, T * result)
    {
        if (v1.size() == v2.size())
        {
            if (!contains<T>(v2, 0))
            {
                #pragma omp simd
                for (int i = 0;i < v1.size();i++)
                {
                    result[i] = (v1[i] / v2[i]);
                }
            }
            else
            {
                throw std::invalid_argument( "Attempted to divide by vector containing zeros");
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector divide");
        }
    }

    //============================================================
    template <class T>
    void add(size_t n, const  T * v, T value, T * result)
    {
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            result[i] = (v[i] + value);
        }
    }

    //============================================================
    template <class T>
    void add(size_t n, const  T * v1, const  T * v2, T * result)
    {
        if (v1.size() == v2.size())
        {            
            int N = v1.size();
            //#pragma omp target map(to:N) map(tofrom: result)
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                result[i] = (v1[i] + v2[i]);
            }         
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector add");
        }
    }

    //============================================================
    template <class T>
    void subtract(size_t n, const  T * v,T value, T * result)
    {        
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            result[i] = (v[i] - value);
        }
    }

    //============================================================
    template <class T>
    void subtract(size_t n, const  T * v1, const  T * v2, T * result)
    {
        if (v1.size() == v2.size())
        {            
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                result[i] = (v1[i] - v2[i]);
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector subtraction");
        }
    }

    //============================================================
    template <class T>
    void abs(size_t n, const  T * v, T * result)
    {        
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            result[i] = std::fabs(v[i]);
        }        
    }

    //============================================================
    template <class T>
    void square(size_t n, const  T * v, T * result)
    {        
        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            result[i] = (v[i]*v[i]);
        }
    }


    //============================================================
    template <class T>
    void squareRoot(size_t n, const  T * v, T * result)
    {
        if (isAllPositiveOrZero(v))
        {            
            #pragma omp simd
            for (int i = 0;i < n;i++)
            {
                result[i] = (T) std::sqrt((T)v[i]);
            }            
        }
        else
        {
            throw std::invalid_argument( "Attempted to compute square root of vector containing negative numbers");
        }
    }

    //============================================================
    template <class T>
    void scale(size_t n, const  T * v,T lowerLimit,T upperLimit, T * result)
    {    
        T minVal = (T) min(v);
        T maxVal = (T) max(v);
        T outputRange = upperLimit - lowerLimit;
        T inputRange = maxVal - minVal;

        #pragma omp simd
        for (int i = 0;i < n;i++)
        {
            T value = (T) v[i];
            T scaledValue = ((value - minVal) * outputRange) / inputRange + lowerLimit;

            result[i] = (scaledValue);
        }
    }

    //============================================================
    template <class T>
    void difference(size_t n, const  T * v, T * result)
    {        
        #pragma omp simd
        for (int i = 1;i < n;i++)
        {
            result[i] = (v[i]-v[i-1]);
        }
    }


    //============================================================
    template <class T>
    void zeros(int N, T * result)
    {
        if (N >= 0)
        {        
            #pragma omp simd
            for (int i = 0;i < N;i++)
            {
                result[i] = T(0);
            }
        }
        else
        {
            throw std::invalid_argument( "Cannot create vector with negative length");
        }
    }

    //============================================================
    template <class T>
    void ones(int N, T * result)
    {
        if (N >= 0)
        {            
            #pragma omp simd
            for (int i = 0;i < N;i++)
            {
                result[i] = T(1);
            }         
        }
        else
        {
            throw std::invalid_argument( "Cannot create vector with negative length");
        }
    }


    //============================================================
    template <class T>
    void range(int limit1,int limit2,int step, T * result)
    {        
        if (step > 0)    
        {
            #pragma omp simd
            for (T i = limit1;i < limit2;i += step)
            {
                result[i] = (i);
            }
        }
        else if (step < 0)
        {
            #pragma omp simd
            for (T i = limit1;i > limit2;i += step)
            {
                result[i] = (i);
            }
        }
        else
        {
            throw std::invalid_argument( "Cannot use a step size of 0 when creating a range of values");
        }
    }

    //============================================================
    template <class T>
    void range(int maxValue, T * result)
    {
        return range<T>(0, maxValue, 1, result);
    }

    //============================================================
    template <class T>
    void range(int minValue,int maxValue, T * result)
    {
        return range<T>(minValue, maxValue, 1, result);
    }

    //============================================================
    template <class T>
    T dotProduct(size_t n, const  T * v1, const  T * v2)
    {
        // if vector size is the same
        if (v1.size() == v2.size())
        {
            T sumVal = 0.0;        
            // sum the element-wise product
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {            
                sumVal += (v1[i]*v2[i]);
            }

            // return the sum as the dot product
            return sumVal;
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in vector dot product");
        }
    }

    //============================================================
    template <class T>
    T euclideanDistance(size_t n, const  T * v1, const  T * v2)
    {
        // if vector size is the same
        if (v1.size() == v2.size())
        {
            T sumVal = 0.0;

            // sum the squared difference
            #pragma omp simd
            for (int i = 0;i < v1.size();i++)
            {
                T diff = (T) (v1[i] - v2[i]);
                sumVal += (diff*diff);
            }

            // if sum is bigger than zero
            if (sumVal > 0)
            {
                // return the square root of the sum as the Euclidean distance
                return std::sqrt(sumVal);
            }
            else // all differences were zero, so report 0.0 as Euclidean distance
            {
                return 0.0;
            }
        }
        else
        {
            throw std::invalid_argument( "Vector lengths differ in Euclidean distance calculation");
        }
    }

    //============================================================
    template <class T>
    T cosineSimilarity(size_t n, const  T * v1, const  T * v2)
    {
    return dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2));
    }

    //============================================================
    template <class T>
    T cosineDistance(size_t n, const  T * v1, const  T * v2)
    {
        return 1.0 - cosineSimilarity(v1, v2);
    }

    
    template<typename T>
    void stdmean(T *sig_src_arr, uint32_t blockSize, T * result){
        T sum = T(0);
        uint32_t blkCnt;
        T in1,in2,in3, in4;
        assert(blockSize != 0);
        //Right shifted by 4 so divided by 4

        #pragma omp simd
        for(blkCnt = blockSize>>2U;blkCnt > 0;blkCnt--){
            in1 = *sig_src_arr++;
            in2 = *sig_src_arr++;
            in3 = *sig_src_arr++;
            in4 = *sig_src_arr++;
            sum += in1;
            sum += in2;
            sum += in3;
            sum += in4;        
        }
        

        #pragma omp simd
        for(blkCnt = blockSize% 0x4;blkCnt > 0;blkCnt--){
            sum += *sig_src_arr++;        
        }
        
        *result = sum/(T)blockSize;        
    }
    

    template<typename T>
    void stdrms(T *pSig_src_arr, uint32_t blockSize, T *pResult)
    {
        T sum = 0.0;
        uint32_t blkCnt;
        T in;
        assert(blockSize != 0);
        
        #pragma omp simd
        for(blkCnt = blockSize >>2;blkCnt > 0;blkCnt--){
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;
            in = *pSig_src_arr++;
            sum += in*in;        
        }
        
        #pragma omp simd
        for(blkCnt = blockSize%0x4;blkCnt>0;blkCnt--)
        {
            in = *pSig_src_arr++;
            sum += in*in;        
        }        
        *pResult = std::sqrt(sum/(T)blockSize);
    }

    
    template<typename T>
    void stddev(T * pSig_src_arr, uint32_t blockSize, T *pResult)
    {

        T sum = 0.0;
        T sumOfSquares = 0.0;
        T in;

        uint32_t blkCnt;

        T meanOfSquares, mean, squareOfMean;
        T squareOfSum = 0.0;

        T var;

        if(blockSize == 1){
            *pResult = 0;
            return;
        }

        

        #pragma omp simd
        for(blkCnt = blockSize>>2;blkCnt>0;blkCnt--){
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;        
        }

        

        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt>0;blkCnt--){
        //perform this operation 4 times
            in = *pSig_src_arr++;
            sum+= in;
            sumOfSquares += in*in;
        }

        meanOfSquares = sumOfSquares / ((T)blockSize-1.0);
        mean = sum/(T) blockSize;

        squareOfMean = (mean*mean) * ((T)blockSize/(T)(blockSize-1.0));

        *pResult = sqrt((meanOfSquares-squareOfMean));
    }

    
    template<typename T>
    void stdvariance(T * pSig_src_arr, uint32_t blockSize, T *pResult)
    {
        T fMean, fValue;
        uint32_t blkCnt;
        T * pInput = pSig_src_arr;

        T sum = 0.0;
        T fSum = 0.0;

        T in1, in2, in3, in4;

        if(blockSize <= 1){
            *pResult = 0;
            return;
        }

        
        #pragma omp simd
        for(blkCnt = blockSize >>2U;blkCnt>0;blkCnt--){
            in1 = *pInput++;
            in2 = *pInput++;
            in3 = *pInput++;
            in4 = *pInput++;

            sum+= in1;
            sum+= in2;
            sum+= in3;
            sum+= in4;    
        }
        

        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt > 0;blkCnt--){
            sum += *pInput++;
        }

        fMean = sum/(T) blockSize;
        pInput = pSig_src_arr;
        
        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt > 0;blkCnt--){
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
        }
        

        #pragma omp simd
        for(blkCnt = blockSize % 0x4;blkCnt>0;blkCnt--){
            fValue = *pInput++ - fMean;
            fSum += fValue*fValue;
        }

        *pResult = fSum/(T)(blockSize-1.0);
    }
}