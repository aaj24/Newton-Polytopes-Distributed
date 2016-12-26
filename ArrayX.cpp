//
//  ArrayX.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 9/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "ArrayX.hpp"
#include <iostream>

template <typename T>
ArrayX<T>::ArrayX(){
    max = 0;
}

template <typename T>
ArrayX<T>::ArrayX(int x){
    max = x;
    array = new T[max];
}

template <typename T>
ArrayX<T>::~ArrayX(){
    delete [] array;
}

template <typename T>
void ArrayX<T>::doubleSize(){
    if (max == 0) {
        max = 1;
        array = new T[max];
    }
    else {
        max <<= 1;
        T* old = array;
        array = new T[max];
        memcpy(array, old, max>>1);
        delete [] old;
    }
}

template <typename T>
void ArrayX<T>::add(T t) {
    if (size + 1 > max)
        doubleSize();
    array[size++] = t;
}

template <typename T>
void ArrayX<T>::empty() {
    size = 0;
}

template <typename T>
T& ArrayX<T>::operator [](int index) {
    return ArrayX<T>::array[index];
}

template <typename T>
ArrayX<T>& ArrayX<T>::operator=(const ArrayX<T>& t) {
    if (this != &t) {
        max = t.max;
        size = t.size;
        delete [] array;
        array = new T[max];
        memcpy(array, t.array, size);
    }
    return *this;
}