# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <mpi.h>
# include <sstream>
# include <string>
# define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()
using namespace std;

int main ( int argc, char *argv[] )
{
	int id;
	int p;
	double wtime;
	MPI::Init ( argc, argv );
	p = MPI::COMM_WORLD.Get_size ( );
	id = MPI::COMM_WORLD.Get_rank ( );
	if ( id == 0 ) 
	{
		cout << "Initiating MPI_BLAST - Master process:\n";
		cout << "  The number of processes is " << p << "\n";
		cout << "\n";
		wtime = MPI::Wtime ( );
	}
	string blast_path = string(argv[1]);
	string allpath = string(argv[2]);
	string cmdStr = blast_path + "blastn -task megablast -db " + allpath + "transHits -evalue 1e-01 -ungapped -num_threads 4  -query " + allpath + "blastQuery_" + SSTR(id) + ".fasta -out " + allpath + "transHitsBLAST_" + SSTR(id) + ".txt -outfmt \"10 qseqid sseqid qstart qend sstart send\" ";
    int ret = system(cmdStr.c_str());
	if ( id == 0 )
	{
		wtime = MPI::Wtime ( ) - wtime;
		cout << "  Elapsed wall clock time = " << wtime << " seconds.\n";
	}
	MPI::Finalize ( );
	return 0;
}
