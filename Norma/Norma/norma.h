//#define N 400


double *W;
double *H;
double *A;

void createArrays(void);
void freeArrays(void);

void solveAnalitic(void);
void solveAnalitic2(void);

void initializeArrays(char *filename);

double normL1();

double normLINF();

void fprintAnalitic();