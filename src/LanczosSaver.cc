#include "LanczosSaver.hh"
bool LanczosSaver::Save(const char * filename,
						const CMatrix * matrix,
						int * r
						)
{
  if(strlen(filename) < 1)
	return false;
  FILE * fout = fopen(filename, "wb");
  if(fout == NULL)
	{
	  throw RLException("Could not open interaction output file.");
	}
  if(matrix == NULL)
	{
	  throw RLException("Cannot save NULL matrix.");
	}

  uint dim = matrix->Rows();
  if(fwrite(&dim, sizeof(dim), 1, fout) != 1)
	throw RLException("File write failed.");

  ulong dimL = (ulong)dim*(ulong)dim;  
  if(fwrite(matrix->GetArray(), sizeof(ComplexDouble), dimL, fout) != dimL)
	throw RLException("File write failed.");

  if(fwrite(&(r[0]), sizeof(int), 2, fout) != 2)
	throw RLException("File write failed.");
  
  ulong validateNo = FILE_VALIDATION_NUMBER_LCS;
  if(fwrite(&validateNo, sizeof(validateNo), 1, fout) != 1)
	throw RLException("File write failed.");

  fclose(fout); 
  fout = NULL;

  return true;
}
