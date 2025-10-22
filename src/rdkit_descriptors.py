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
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from descriptor_generator import AbstractCalculator, get_base_parser


def main():
    # Usage:
    # python src/rdkit_descriptors.py -i data/10.smi -o descriptors-rdkit.smi -d tab --write-header

    # ----- command line args definitions ---------------------------------------------

    parser = get_base_parser()

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    descriptors = [
        d[0] for d in Descriptors._descList  # pylint: disable=protected-access
    ]
    desc_calculator = MoleculeDescriptors.MolecularDescriptorCalculator(
        descriptors,
        missing_val=args.missing_val,
    )

    class DescriptorCalculator(AbstractCalculator):
        calculator = desc_calculator.CalcDescriptors
        descriptor_names = descriptors

    calc = DescriptorCalculator(args)
    print("got calc", calc)
    calc.run()


if __name__ == "__main__":
    main()
