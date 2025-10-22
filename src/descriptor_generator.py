# Copyright 2024 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Generates rdkit descriptors for molecules.

import sys
import time
from abc import ABC
from typing import Callable

from dm_job_utilities.dm_log import DmLog
from rdkit import Chem

import rdkit_utils
import utils


class AbstractCalculator(ABC):

    calculator: Callable
    descriptor_names: list[str]

    def __init__(self, arguments):
        self.filename = arguments.input
        self.outfile = arguments.output
        self.mode = arguments.fragment_method
        self.include_3d = arguments.include_3d
        self.omit_fields = arguments.omit_fields
        self.id_column = arguments.id_column
        self.mol_column = arguments.mol_column
        self.read_header = arguments.read_header
        self.write_header = arguments.write_header
        self.read_records = arguments.read_records
        self.interval = arguments.interval
        self.missing_val = arguments.missing_val
        self.delimiter = utils.read_delimiter(arguments.delimiter)

    def __init_subclass__(cls):
        super().__init_subclass__()
        if not hasattr(cls, "calculator"):
            raise TypeError(f"{cls.__name__} must define 'calculator'")
        if not hasattr(cls, "descriptor_names"):
            raise TypeError(f"{cls.__name__} must define 'descriptor_names'")

    def run(
        self,
        # filename,
        # outfile,
        # mode="hac",
        # include_3d=False,
        # delimiter=None,
        # id_column=None,
        # mol_column=0,
        # omit_fields=False,
        # read_header=False,
        # write_header=False,
        # read_records=50,
        # interval=1000,
        # missing_val=None,
    ):

        if self.include_3d:
            if not (
                self.filename.endswith(".sdf") or self.filename.endswith(".sdf.gz")
            ):
                DmLog.emit_event(
                    "If calculating 3D descriptors input must be a SDF with 3D molecules"
                )
                exit(1)
            if not rdkit_utils.check_molecules_are_3d(self.filename):
                DmLog.emit_event("Molecules do not seem to have 3D coordinates")
                exit(1)

        utils.expand_path(self.outfile)

        count = 0
        errors = 0

        t0 = time.time()

        # setup the reader
        reader = rdkit_utils.create_reader(
            self.filename,
            id_column=self.id_column,
            mol_column=self.mol_column,
            read_records=self.read_records,
            read_header=self.read_header,
            delimiter=self.delimiter,
        )
        extra_field_names = reader.get_extra_field_names()

        # calc_field_names = [
        #     d[0] for d in Descriptors._descList  # pylint: disable=protected-access
        # ]
        # calc = MoleculeDescriptors.MolecularDescriptorCalculator(calc_field_names)

        # # properties in header need to be in same order as later calculated
        # calc_morganfp = None
        # if morgan_fingerprints:
        #     # redundant?
        #     # morganfp_as_bitvect = morgan_fingerprints.pop('as_bitvect')
        #     # if morganfp_as_bitvect:
        #     #     calc_field_names.append('RDKFingerprintBitVect')
        #     # else:
        #     #     calc_field_names.append('RDKFingerprint')
        #     calc_field_names.append("MorganFingerprint")
        #     calc_morganfp = rdFingerprintGenerator.GetMorganGenerator(**morgan_fingerprints)

        # calc_rdkitfp = None
        # if rdkit_fingerprints:
        #     calc_rdkitfp = rdFingerprintGenerator.GetRDKitFPGenerator(**rdkit_fingerprints)
        #     calc_field_names.append("RDKFingerprint")

        DmLog.emit_event(
            "Calculating {} descriptors".format(len(self.descriptor_names))
        )

        # setup the writer
        writer = rdkit_utils.create_writer(
            self.outfile,
            extra_field_names=extra_field_names,
            calc_prop_names=self.descriptor_names,
            delimiter=self.delimiter,
            id_column=self.id_column,
            mol_column=self.mol_column,
        )

        id_col_type, id_col_value = utils.is_type(self.id_column, int)

        # read the input records and write the output
        while True:
            t = reader.read()
            # break if no more data to read
            if not t:
                break
            mol, smi, mol_id, props = t

            if count == 0 and self.write_header:
                headers = rdkit_utils.generate_headers(
                    id_col_type,
                    id_col_value,
                    reader.get_mol_field_name(),
                    reader.field_names,
                    self.descriptor_names,
                    self.omit_fields,
                )

                writer.write_header(headers)

            count += 1

            if self.interval and count % self.interval == 0:
                DmLog.emit_event("Processed {} records".format(count))
            if count % 50000 == 0:
                # Emit a 'total' cost, replacing all prior costs
                DmLog.emit_cost(count)

            if not mol:
                errors += 1
                DmLog.emit_event("Failed to read molecule for record", count)
                continue

            if self.omit_fields:
                for name in mol.GetPropNames():
                    mol.ClearProp(name)
                props = []

            try:
                if self.mode == "none":
                    biggest = mol
                    cann_smi = smi
                else:
                    biggest = rdkit_utils.fragment(mol, self.mode)
                    cann_smi = Chem.MolToSmiles(biggest)

                # this is now the calculator defined in child classes
                values = self.calculator(biggest)

            except KeyboardInterrupt:
                utils.log("Interrupted")
                sys.exit(0)

            except Exception as exc:
                print(exc)
                errors += 1
                DmLog.emit_event("Failed to process record", count)
                continue

            writer.write(cann_smi, biggest, mol_id, props, values)

        writer.close()
        reader.close()

        t1 = time.time()
        utils.log("Processing took {} secs".format(round(t1 - t0)))


# def run(
#     filename,
#     outfile,
#     mode="hac",
#     include_3d=False,
#     delimiter=None,
#     id_column=None,
#     mol_column=0,
#     omit_fields=False,
#     read_header=False,
#     write_header=False,
#     read_records=50,
#     interval=1000,
#     missing_val=None,
#     rdkit_fingerprints=None,
#     morgan_fingerprints=None,
# ):

#     if include_3d:
#         if not (filename.endswith(".sdf") or filename.endswith(".sdf.gz")):
#             DmLog.emit_event(
#                 "If calculating 3D descriptors input must be a SDF with 3D molecules"
#             )
#             exit(1)
#         if not rdkit_utils.check_molecules_are_3d(filename):
#             DmLog.emit_event("Molecules do not seem to have 3D coordinates")
#             exit(1)

#     utils.expand_path(outfile)

#     count = 0
#     errors = 0

#     t0 = time.time()

#     # setup the reader
#     reader = rdkit_utils.create_reader(
#         filename,
#         id_column=id_column,
#         mol_column=mol_column,
#         read_records=read_records,
#         read_header=read_header,
#         delimiter=delimiter,
#     )
#     extra_field_names = reader.get_extra_field_names()

#     calc_field_names = [
#         d[0] for d in Descriptors._descList  # pylint: disable=protected-access
#     ]
#     calc = MoleculeDescriptors.MolecularDescriptorCalculator(calc_field_names)

#     # properties in header need to be in same order as later calculated
#     calc_morganfp = None
#     if morgan_fingerprints:
#         # redundant?
#         # morganfp_as_bitvect = morgan_fingerprints.pop('as_bitvect')
#         # if morganfp_as_bitvect:
#         #     calc_field_names.append('RDKFingerprintBitVect')
#         # else:
#         #     calc_field_names.append('RDKFingerprint')
#         calc_field_names.append("MorganFingerprint")
#         calc_morganfp = rdFingerprintGenerator.GetMorganGenerator(**morgan_fingerprints)

#     calc_rdkitfp = None
#     if rdkit_fingerprints:
#         calc_rdkitfp = rdFingerprintGenerator.GetRDKitFPGenerator(**rdkit_fingerprints)
#         calc_field_names.append("RDKFingerprint")

#     DmLog.emit_event("Calculating {} descriptors".format(len(calc_field_names)))

#     # setup the writer
#     writer = rdkit_utils.create_writer(
#         outfile,
#         extra_field_names=extra_field_names,
#         calc_prop_names=calc_field_names,
#         delimiter=delimiter,
#         id_column=id_column,
#         mol_column=mol_column,
#     )

#     id_col_type, id_col_value = utils.is_type(id_column, int)

#     # read the input records and write the output
#     while True:
#         t = reader.read()
#         # break if no more data to read
#         if not t:
#             break
#         mol, smi, mol_id, props = t

#         if count == 0 and write_header:
#             headers = rdkit_utils.generate_headers(
#                 id_col_type,
#                 id_col_value,
#                 reader.get_mol_field_name(),
#                 reader.field_names,
#                 calc_field_names,
#                 omit_fields,
#             )

#             writer.write_header(headers)

#         count += 1

#         if interval and count % interval == 0:
#             DmLog.emit_event("Processed {} records".format(count))
#         if count % 50000 == 0:
#             # Emit a 'total' cost, replacing all prior costs
#             DmLog.emit_cost(count)

#         if not mol:
#             errors += 1
#             DmLog.emit_event("Failed to read molecule for record", count)
#             continue

#         if omit_fields:
#             for name in mol.GetPropNames():
#                 mol.ClearProp(name)
#             props = []

#         try:
#             if mode == "none":
#                 biggest = mol
#                 cann_smi = smi
#             else:
#                 biggest = rdkit_utils.fragment(mol, mode)
#                 cann_smi = Chem.MolToSmiles(biggest)

#             values = calc.CalcDescriptors(biggest, missing_val=missing_val)
#             # print(values)

#             # preselect? but signatures are different..
#             # should there be different options in the same output? like nBits = [1024, 2048,..]
#             # print(morgan_fingerprints)
#             values = list(values)
#             if calc_morganfp:
#                 fp = calc_morganfp.GetFingerprint(mol)
#                 values.append(fp.ToBitString())

#             if calc_rdkitfp:
#                 fp = calc_rdkitfp.GetFingerprint(mol)
#                 values.append(fp.ToBitString())

#             # if True:
#             #     fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
#             #     print(fp)

#         except KeyboardInterrupt:
#             utils.log("Interrupted")
#             sys.exit(0)

#         except Exception as exc:
#             print(exc)
#             errors += 1
#             DmLog.emit_event("Failed to process record", count)
#             continue

#         writer.write(cann_smi, biggest, mol_id, props, values)

#     writer.close()
#     reader.close()

#     t1 = time.time()
#     utils.log("Processing took {} secs".format(round(t1 - t0)))


# def _obsolete_main():

#     # Run using conda env created from environment-im-mordred.yaml
#     #   or docker environment created from Dockerfile-mordred

#     # Examples:
#     #   python -m im_mordred.descriptor_generator -i data/10.smi -o descriptors.smi -d tab
#     #   python -m im_mordred.descriptor_generator -i data/10+H.smi -o descriptors.smi -d tab --id_column 1 --read-header --writeheader
#     # python src/descriptor_generator.py -i data/10.smi -o descriptors.smi -d tab --write-header --rdkfp-calc-rdkfp --morganfp-calc-morganfp

#     # ----- command line args definitions ---------------------------------------------

#     parser = argparse.ArgumentParser(description="Mordred 2D descriptors")
#     input_group = parser.add_argument_group("Input/output options")
#     input_group.add_argument(
#         "-i", "--input", required=True, help="Input file (.smi or .sdf)"
#     )
#     input_group.add_argument(
#         "-o", "--output", default="descriptors2d.smi", help="Output file (.smi or .sdf"
#     )
#     input_group.add_argument(
#         "--omit-fields",
#         action="store_true",
#         help="Don't include fields from the input in the output",
#     )

#     # to pass tab as the delimiter specify it as $'\t' or use one of
#     # the symbolic names 'comma', 'tab', 'space' or 'pipe'
#     input_group.add_argument("-d", "--delimiter", help="Delimiter when using SMILES")
#     input_group.add_argument(
#         "--id-column",
#         help="Column for name field (zero based integer for .smi, text for SDF)",
#     )
#     input_group.add_argument(
#         "--mol-column",
#         type=int,
#         default=0,
#         help="Column index for molecule when using delineated text formats (zero based integer)",
#     )
#     input_group.add_argument(
#         "--read-header",
#         action="store_true",
#         help="Read a header line with the field names when reading .smi or .txt",
#     )
#     input_group.add_argument(
#         "--write-header",
#         action="store_true",
#         help="Write a header line when writing .smi or .txt",
#     )
#     input_group.add_argument(
#         "--read-records",
#         default=100,
#         type=int,
#         help="Read this many records to determine the fields that are present",
#     )
#     input_group.add_argument(
#         "--interval", default=1000, type=int, help="Reporting interval"
#     )

#     rdkit_generic_group = parser.add_argument_group("General RDKit options")
#     rdkit_generic_group.add_argument(
#         "--missing-val",
#         default=None,
#         help="Missing value",
#     )
#     rdkit_generic_group.add_argument(
#         "--fragment-method",
#         choices=["hac", "mw", "none"],
#         default="hac",
#         help="Strategy for picking largest fragment (mw or hac or none",
#     )

#     rdkit_generic_group.add_argument(
#         "--include-3d",
#         action="store_true",
#         help="Include 3D descriptors (requires 3D molecules in SDF file)",
#     )

#     rdkfp_group = parser.add_argument_group("RDKFingerprint options")
#     rdkfp_group.add_argument(
#         "--rdkfp-calc-rdkfp", action="store_true", help="Calculate RDKFingerprint"
#     )
#     rdkfp_group.add_argument(
#         "--rdkfp-min-path",
#         default=1,
#         type=int,
#         help="Minimum number of bonds to include in the subgraphs",
#     )
#     rdkfp_group.add_argument(
#         "--rdkfp-max-path",
#         default=7,
#         type=int,
#         help="Maximum number of bonds to include in the subgraphs",
#     )
#     rdkfp_group.add_argument(
#         "--rdkfp-size",
#         default=2048,
#         type=int,
#         help="Number of bits in the fingerprint",
#     )
#     # rdkfp_group.add_argument(
#     #     "--rdkfp-bits-per-hash",
#     #     default=2,
#     #     type=int,
#     #     help="Number of bits to set per path",
#     # )
#     rdkfp_group.add_argument(
#         "--rdkfp-bits-per-feature",
#         default=2,
#         type=int,
#         help="Number of bits to set per path",
#     )
#     rdkfp_group.add_argument(
#         "--rdkfp-use-hydrogens",
#         # default=True,
#         # type=bool,
#         action="store_true",
#         help="Include paths involving Hs in the fingerprint if the molecule has explicit Hs",
#     )
#     # rdkfp_group.add_argument(
#     #     "--rdkfp-tgt-density",
#     #     default=0.0,
#     #     type=float,
#     #     help="Fold the fingerprint until this minimum density has been reached",
#     # )
#     # rdkfp_group.add_argument(
#     #     "--rdkfp-min-size",
#     #     default=128,
#     #     type=int,
#     #     help="The minimum size the fingerprint will be folded to when trying to reach tgtDensity",
#     # )
#     rdkfp_group.add_argument(
#         "--rdkfp-branched-paths",
#         # default=True,
#         # type=bool,
#         action="store_true",
#         help="If set, both branched and unbranched paths will be used in the fingerprint",
#     )
#     rdkfp_group.add_argument(
#         "--rdkfp-use-bond-order",
#         # default=True,
#         # type=bool,
#         action="store_true",
#         help="If set, both bond orders will be used in the path hashes",
#     )
#     rdkfp_group.add_argument(
#         "--rdkfp-count-simulation",
#         # default=False,
#         # type=bool,
#         action="store_true",
#         help="If set, both bond orders will be used in the path hashes",
#     )
#     # I don't think it's practical to expose atom pair parameters

#     # where help is commented, cannot find docs
#     morganfp_group = parser.add_argument_group("Morgan fingerprint options")
#     morganfp_group.add_argument(
#         "--morganfp-calc-morganfp",
#         action="store_true",
#         help="Calculate Morgan Fingerprint",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-radius",
#         default=3,
#         type=int,
#         help="Minimum number of bonds to include in the subgraphs",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-count-simulation",
#         # default=False,
#         # type=bool,
#         action="store_true",
#         help="If set, use count simulation while generating the fingerprint",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-use-chirality",
#         # default=False,
#         # type=bool,
#         action="store_true",
#         help="If set, chirality information will be added to the generated fingerprint",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-use-bond-types",
#         # default=True,
#         # type=bool,
#         action="store_true",
#         help="If set, bond types will be included as a part of the default bond invariants",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-only-nonzero-invariants",
#         # default=False,
#         # type=bool,
#         action="store_true",
#         # help="If set, bond types will be included as a part of the default bond invariants",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-include-ring-membership",
#         # default=True,
#         # type=bool,
#         action="store_true",
#         # help="If set, bond types will be included as a part of the default bond invariants",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-count-bonds",
#         # default=False,
#         # type=bool,
#         action="store_true",
#         help="Boundaries for count simulation, corresponding bit will be set if the count is higher than the number provided for that spot",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-size",
#         default=2048,
#         type=int,
#         help="Size of the generated fingerprint, does not affect the sparse versions",
#     )
#     morganfp_group.add_argument(
#         "--morganfp-include-redundandt-environments",
#         # default=False,
#         # type=bool,
#         action="store_true",
#         help="Include redundant environments in the fingerprint",
#     )
#     # morganfp_group.add_argument(
#     #     "--morganfp-as-bitvect",
#     #     default=False,
#     #     type=bool,
#     #     # help="Minimum number of bonds to include in the subgraphs",
#     # )

#     # TODO:
#     # apgen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
#     # ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)

#     # rdkit_descriptor_group = parser.add_argument_group("RDKit descriptors options")

#     args = parser.parse_args()
#     DmLog.emit_event("descriptor_calc: ", args)

#     delimiter = utils.read_delimiter(args.delimiter)

#     if args.rdkfp_calc_rdkfp:
#         rdkit_fingerprints = {
#             "minPath": args.rdkfp_min_path,
#             "maxPath": args.rdkfp_max_path,
#             "fpSize": args.rdkfp_size,
#             # "nBitsPerHash": args.rdkfp_bits_per_hash,
#             "numBitsPerFeature": args.rdkfp_bits_per_feature,
#             "useHs": args.rdkfp_use_hydrogens,
#             # "tgtDensity": args.rdkfp_tgt_density,
#             # "minSize": args.rdkfp_min_size,
#             "branchedPaths": args.rdkfp_branched_paths,
#             "useBondOrder": args.rdkfp_use_bond_order,
#             "countSimulation": args.rdkfp_count_simulation,
#         }
#     else:
#         rdkit_fingerprints = {}

#     if args.morganfp_calc_morganfp:
#         morgan_fingerprints = {
#             "radius": args.morganfp_radius,
#             "countSimulation": args.morganfp_count_simulation,
#             "includeChirality": args.morganfp_use_chirality,
#             "useBondTypes": args.morganfp_use_bond_types,
#             "onlyNonzeroInvariants": args.morganfp_only_nonzero_invariants,
#             "includeRingMembership": args.morganfp_include_ring_membership,
#             "fpSize": args.morganfp_size,
#             "includeRedundantEnvironments": args.morganfp_include_redundandt_environments,
#             # "as_bitvect": args.morganfp_as_bitvect,
#         }
#     else:
#         morgan_fingerprints = {}

#     # run(
#     #     args.input,
#     #     args.output,
#     #     mode=args.fragment_method,
#     #     include_3d=args.include_3d,
#     #     omit_fields=args.omit_fields,
#     #     delimiter=delimiter,
#     #     id_column=args.id_column,
#     #     mol_column=args.mol_column,
#     #     read_header=args.read_header,
#     #     write_header=args.write_header,
#     #     read_records=args.read_records,
#     #     interval=args.interval,
#     #     missing_val=args.missing_val,
#     #     rdkit_fingerprints=rdkit_fingerprints,
#     #     morgan_fingerprints=morgan_fingerprints,
#     # )


# # if __name__ == "__main__":
# #     main()


# def get_base_parser():

#     # Run using conda env created from environment-im-mordred.yaml
#     #   or docker environment created from Dockerfile-mordred

#     # Examples:
#     #   python -m im_mordred.descriptor_generator -i data/10.smi -o descriptors.smi -d tab
#     #   python -m im_mordred.descriptor_generator -i data/10+H.smi -o descriptors.smi -d tab --id_column 1 --read-header --writeheader
#     # python src/descriptor_generator.py -i data/10.smi -o descriptors.smi -d tab --write-header --rdkfp-calc-rdkfp --morganfp-calc-morganfp

#     # ----- command line args definitions ---------------------------------------------

#     parser = argparse.ArgumentParser(description="Mordred 2D descriptors")
#     input_group = parser.add_argument_group("Input/output options")
#     input_group.add_argument(
#         "-i", "--input", required=True, help="Input file (.smi or .sdf)"
#     )
#     input_group.add_argument(
#         "-o", "--output", default="descriptors2d.smi", help="Output file (.smi or .sdf"
#     )
#     input_group.add_argument(
#         "--omit-fields",
#         action="store_true",
#         help="Don't include fields from the input in the output",
#     )

#     # to pass tab as the delimiter specify it as $'\t' or use one of
#     # the symbolic names 'comma', 'tab', 'space' or 'pipe'
#     input_group.add_argument("-d", "--delimiter", help="Delimiter when using SMILES")
#     input_group.add_argument(
#         "--id-column",
#         help="Column for name field (zero based integer for .smi, text for SDF)",
#     )
#     input_group.add_argument(
#         "--mol-column",
#         type=int,
#         default=0,
#         help="Column index for molecule when using delineated text formats (zero based integer)",
#     )
#     input_group.add_argument(
#         "--read-header",
#         action="store_true",
#         help="Read a header line with the field names when reading .smi or .txt",
#     )
#     input_group.add_argument(
#         "--write-header",
#         action="store_true",
#         help="Write a header line when writing .smi or .txt",
#     )
#     input_group.add_argument(
#         "--read-records",
#         default=100,
#         type=int,
#         help="Read this many records to determine the fields that are present",
#     )
#     input_group.add_argument(
#         "--interval", default=1000, type=int, help="Reporting interval"
#     )

#     rdkit_generic_group = parser.add_argument_group("General RDKit options")
#     rdkit_generic_group.add_argument(
#         "--missing-val",
#         default=None,
#         help="Missing value",
#     )
#     rdkit_generic_group.add_argument(
#         "--fragment-method",
#         choices=["hac", "mw", "none"],
#         default="hac",
#         help="Strategy for picking largest fragment (mw or hac or none",
#     )

#     rdkit_generic_group.add_argument(
#         "--include-3d",
#         action="store_true",
#         help="Include 3D descriptors (requires 3D molecules in SDF file)",
#     )


#     return parser
#     # delimiter = utils.read_delimiter(args.delimiter)


#     # run(
#     #     args.input,
#     #     args.output,
#     #     mode=args.fragment_method,
#     #     include_3d=args.include_3d,
#     #     omit_fields=args.omit_fields,
#     #     delimiter=delimiter,
#     #     id_column=args.id_column,
#     #     mol_column=args.mol_column,
#     #     read_header=args.read_header,
#     #     write_header=args.write_header,
#     #     read_records=args.read_records,
#     #     interval=args.interval,
#     #     missing_val=args.missing_val,
#     #     rdkit_fingerprints=rdkit_fingerprints,
#     #     morgan_fingerprints=morgan_fingerprints,
#     # )
