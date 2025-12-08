#!/usr/bin/env python

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


from dm_job_utilities.dm_log import DmLog
from rdkit.Chem import rdFingerprintGenerator

from descriptor_generator import AbstractCalculator, get_base_parser


def main():

    # Usage:
    # python src/fp_rdkit.py -i data/10.smi -o descriptors-rdkfp.smi -d tab --write-header

    # ----- command line args definitions ---------------------------------------------

    parser = get_base_parser()

    rdkfp_group = parser.add_argument_group("RDKFingerprint options")
    # rdkfp_group.add_argument(
    #     "--rdkfp-calc-rdkfp", action="store_true", help="Calculate RDKFingerprint"
    # )
    rdkfp_group.add_argument(
        "--rdkfp-min-path",
        default=1,
        type=int,
        help="Minimum number of bonds to include in the subgraphs",
    )
    rdkfp_group.add_argument(
        "--rdkfp-max-path",
        default=7,
        type=int,
        help="Maximum number of bonds to include in the subgraphs",
    )
    rdkfp_group.add_argument(
        "--rdkfp-size",
        default=2048,
        type=int,
        help="Number of bits in the fingerprint",
    )
    # rdkfp_group.add_argument(
    #     "--rdkfp-bits-per-hash",
    #     default=2,
    #     type=int,
    #     help="Number of bits to set per path",
    # )
    rdkfp_group.add_argument(
        "--rdkfp-bits-per-feature",
        default=2,
        type=int,
        help="Number of bits to set per path",
    )
    rdkfp_group.add_argument(
        "--rdkfp-use-hydrogens",
        # default=True,
        # type=bool,
        action="store_true",
        help="Include paths involving Hs in the fingerprint if the molecule has explicit Hs",
    )
    # rdkfp_group.add_argument(
    #     "--rdkfp-tgt-density",
    #     default=0.0,
    #     type=float,
    #     help="Fold the fingerprint until this minimum density has been reached",
    # )
    # rdkfp_group.add_argument(
    #     "--rdkfp-min-size",
    #     default=128,
    #     type=int,
    #     help="The minimum size the fingerprint will be folded to when trying to reach tgtDensity",
    # )
    rdkfp_group.add_argument(
        "--rdkfp-branched-paths",
        # default=True,
        # type=bool,
        action="store_true",
        help="If set, both branched and unbranched paths will be used in the fingerprint",
    )
    rdkfp_group.add_argument(
        "--rdkfp-use-bond-order",
        # default=True,
        # type=bool,
        action="store_true",
        help="If set, both bond orders will be used in the path hashes",
    )
    rdkfp_group.add_argument(
        "--rdkfp-count-simulation",
        # default=False,
        # type=bool,
        action="store_true",
        help="If set, both bond orders will be used in the path hashes",
    )
    # I don't think it's practical to expose atom pair parameters

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    rdkit_fingerprints = {
        "minPath": args.rdkfp_min_path,
        "maxPath": args.rdkfp_max_path,
        "fpSize": args.rdkfp_size,
        # "nBitsPerHash": args.rdkfp_bits_per_hash,
        "numBitsPerFeature": args.rdkfp_bits_per_feature,
        "useHs": args.rdkfp_use_hydrogens,
        # "tgtDensity": args.rdkfp_tgt_density,
        # "minSize": args.rdkfp_min_size,
        "branchedPaths": args.rdkfp_branched_paths,
        "useBondOrder": args.rdkfp_use_bond_order,
        "countSimulation": args.rdkfp_count_simulation,
    }

    fpcalculator = rdFingerprintGenerator.GetRDKitFPGenerator(**rdkit_fingerprints)

    class RDFPCalculator(AbstractCalculator):
        calculator = fpcalculator.GetFingerprint
        descriptor_names = [
            "RDKFingerprint",
        ]

        def calculate(self, *args):
            result = self.calculator(*args)
            return (result.ToBitString(),)

    calc = RDFPCalculator(args)
    calc.run()


if __name__ == "__main__":
    main()
