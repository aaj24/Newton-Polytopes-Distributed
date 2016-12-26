//
//  ArrayX.hpp
//  CMPS 299
//
//  Created by Ali Jbaily on 9/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#ifndef ArrayX_hpp
#define ArrayX_hpp

#include <utility>
#include <string.h>

template <typename T>
class ArrayX {
    T* array;
    int max;
    
public:
    int size = 0;
    void add(T);
    void empty();
    ArrayX<T>();
    ArrayX<T>(int);
    ~ArrayX();
    ArrayX(const ArrayX &obj);  // copy constructor
    T& operator [](int);
    ArrayX& operator=(const ArrayX&);
    
private:
    void doubleSize();
};

#endif /* ArrayX_hpp */
