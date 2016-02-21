// Test function
void hello();
//////////////////////////////////////////////////////////
// histos
int decla(int,int);
int save_hist(int);
void style();
void draw();
void draw_out();
void load();
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// tags
double GenLevelWeight(vector<PseudoJet>, int,int);
int HH4b(vector<PseudoJet> jets, vector<int> btag, vector<int> btrue, double weight,int cluster, int);
////
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets, vector<int> & btag, 
             vector<int> & bmistag, vector<int> & fattag, vector<int> & btrue, double weight, int cluster, int);
int recol(vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, double weight);
bool recohadt(int & bh, int & bl, vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type);
bool recolept2step(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
                   vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type);
bool recotlepeq(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
                vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type);

int fullylep(int & bh, int & bl, int & ,  vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met,double weight);
