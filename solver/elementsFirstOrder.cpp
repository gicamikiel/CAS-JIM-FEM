#include <elementsFirstOrder.hpp>
#include <stdexcept>

int FirstOrderQuad::getNumNodes() const {
    return 4;
}

double FirstOrderQuad::getXi(int node) const {
    return parXi[node];
}

double FirstOrderQuad::getEta(int node) const {
    return parEta[node];
}

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
    throw std::invalid_argument("invalid shape function subscript (local node number).");
}

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
    throw std::invalid_argument("invalid shape function subscript (local node number) or coordinate.");
}