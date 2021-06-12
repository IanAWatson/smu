// Generate small molecules exhaustively

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/molecule.h"

#include "smu.h"

namespace smu {

int verbose = 0;

Report_Progress report_progress;

void
usage(int rc) {
  cerr << "Generates SMU molecules\n";
  cerr << " -M <number>       max number of atoms in generated molecules\n";
  cerr << " -x                use non aromatic unique smiles\n";
  cerr << " -O                allow O- to be added to any atom (not just N+)\n";
  cerr << " -f <sep>          output separator\n";
  exit(0);
}

int generate_smu(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vA:E:M:xOR:f:");

  if (cl.unrecognised_options_encountered())
  {
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! process_elements(cl, verbose))
  {
    cerr << "Cannot parse element specifications\n";
    usage(2);
  }

  if (cl.option_present('R')) {
    if (! report_progress.initialise(cl, 'R', verbose)) {
      cerr << "Cannot initialise progress reporting (-R)\n";
      return 1;
    }
  }

  SMU_Params params;
  params.to_add.add(new Atom(6));
  params.to_add.add(new Atom(7));
  params.to_add.add(new Atom(7));
  params.to_add.last_item()->set_formal_charge(1);
  params.to_add.add(new Atom(8));
  params.to_add.add(new Atom(8));
  params.to_add.last_item()->set_formal_charge(-1);
  params.to_add.add(new Atom(9));

  int max_atoms = 7;
  if (cl.option_present('M')) {
    if (! cl.value('M', max_atoms) || max_atoms < 2) {
      cerr << "Invalid max atoms (-M)\n";
      usage(1);
    }
  }

  params.max_atoms = max_atoms;

  params.oneg_only_added_to_npos = true;
  if (cl.option_present('O')) {
    params.oneg_only_added_to_npos = false;
    if (verbose)
      cerr << "O- can be added to any existing atom\n";
  }

  params.non_aromatic_unique_smiles = false;
  if (cl.option_present('x')) {
    params.non_aromatic_unique_smiles = true;
    if (verbose)
      cerr << "Will use non aromatic unique smiles\n";
  }

  Molecule carbon;
  carbon.build_from_smiles("C");

  SmuResults smu_results;

  if (cl.option_present('f')) {
    IWString f = cl.option_value('f');
    if (f == "space") {
      f = ' ';
    } else if (f == "tab") {
      f = '\t';
    }
    smu_results.set_output_separator(f[0]);
  }

  Expand(carbon, params, smu_results);
  if (verbose) {
    smu_results.Report(std::cerr);
    cerr << "Generated " << smu_results.size() << " molecules\n";
  }

  smu_results.Write(std::cout);

  return 0;
}

}  // namespace smu

int
main(int argc, char ** argv) {
  return smu::generate_smu(argc, argv);
}
