#ifndef EVOLUTIONMODEL_H_INCLUDED
#define EVOLUTIONMODEL_H_INCLUDED

class EvolutionModel {
  public:
    virtual double P(char a, char b, double t) const = 0;
    virtual double prior(char s) const = 0;
};

class Kimura : public EvolutionModel {
public:
    Kimura(double R) { this->R = R; }
    virtual double P(char a, char b, double t) const;
    double prior(char s) const;
    void setR(double R) { this->R = R; }
private:
    double R;
};

class JukesCantor : public EvolutionModel {
public:
    JukesCantor(double u) { this->u = u; }
    virtual double P(char a, char b, double t) const;
    double prior(char s) const;
    void setU(double u) { this->u = u; }
private:
    double u;
};

#endif
