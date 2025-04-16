#include <elementsFirstOrder.hpp>
#include <stdexcept>

// FirstOrderQuad
/*
For a quad, assume the following node arrangement
3--2  ^ eta
|  |  | 
0--1  --> xi
*/

KOKKOS_FUNCTION
int FirstOrderQuad::getNumNodes() const {
    return 4;
}

KOKKOS_FUNCTION
double FirstOrderQuad::getXi(int node) const {
    return parXi[node];
} 

KOKKOS_FUNCTION
double FirstOrderQuad::getEta(int node) const {
    return parEta[node];
}

KOKKOS_FUNCTION
double FirstOrderQuad::shape(int node, double xi, double eta) const {
    switch(node) {
        case 0:
            return 0.25*(1-xi)*(1-eta);
        case 1:
            return 0.25*(1+xi)*(1-eta);
        case 2: 
            return 0.25*(1+xi)*(1+eta);
        case 3:
            return 0.25*(1-xi)*(1+eta);
    }
    return 0;
}

KOKKOS_FUNCTION
double FirstOrderQuad::shapeD(int node, double xi, double eta, CoordinateParametric dim) const {
    switch(node) {
        case 0:
            switch(dim) {
                case XI:
                    return 0.25*-1*(1-eta);
                case ETA:
                    return 0.25*(1-xi)*-1;
            }
        case 1:
            switch(dim) {
                case XI:
                    return 0.25*1*(1-eta);
                case ETA:
                    return 0.25*(1+xi)*-1;
            }
        case 2: 
            switch(dim) {
                case XI:
                    return 0.25*1*(1+eta);
                case ETA:
                    return 0.25*(1+xi)*+1;
            }
        case 3:
            switch(dim) {
                case XI:
                    return 0.25*-1*(1+eta);
                case ETA:
                    return 0.25*(1-xi)*+1;
            }
    }
    return 0;
}

// FirstOrderTri
// For a tri, assume the following node arrangement
/*
3
| \     ^ eta
|   \   |
0----2  --> xi
*/

KOKKOS_FUNCTION
int FirstOrderTri::getNumNodes() const {
    return 3;
}

KOKKOS_FUNCTION
double FirstOrderTri::getXi(int node) const {
    return parXi[node];
}

double FirstOrderTri::getEta(int node) const {
    return parEta[node];
}

KOKKOS_FUNCTION
double FirstOrderTri::shape(int node, double xi, double eta) const {
    switch(node) {
        case 0:
            return 1-xi-eta;
        case 1:
            return xi;
        case 2: 
            return eta;
    }
    return 0;
}

KOKKOS_FUNCTION
double FirstOrderTri::shapeD(int node, double xi, double eta, CoordinateParametric dim) const {
    switch(node) {
        case 0:
            return -1;
        case 1:
            switch(dim) {
                case XI:
                    return 1;
                case ETA:
                    return 0;
            }
        case 2: 
            switch(dim) {
                case XI:
                    return 0;
                case ETA:
                    return 1;
            }
    }
    return 0;
}