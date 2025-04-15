#ifndef SOLVER_ELEMENT
#define SOLVER_ELEMENT

enum CoordinateParametric {
    NONE = 0,
    XI = 1,
    ETA = 2
};

class FiniteElement {
    public:
        virtual int getNumNodes() const;
        virtual double getXi(int idx) const;
        virtual double getEta(int idx) const;
        virtual double shape(int idx, double xi, double eta) const;
        virtual double shapeD(int idx, double xi, double eta, CoordinateParametric dim) const;

    private:
};

#endif