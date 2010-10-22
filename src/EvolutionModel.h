#ifndef EVOLUTIONMODEL_H_INCLUDED
#define EVOLUTIONMODEL_H_INCLUDED

class EvolutionModel {
  public:
    virtual double P(char a, char b, double t) const = 0;
};

class Kimura : public EvolutionModel {
public:
    Kimura(double R) { this->R = R; }
    virtual double P(char a, char b, double t) const;
    void setR(double R) { this->R = R; }
private:
    double R;
};

#endif
