-- Copyright [1999-2014] EMBL-European Bioinformatics Institute
-- and Wellcome Trust Sanger Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

-- Remove Submit* analysis entries

DELETE FROM analysis WHERE logic_name like 'Submit%';

-- Remove pipeline tables

DROP TABLE input_id_analysis;

DROP TABLE input_id_type_analysis;

DROP TABLE job;

DROP TABLE job_status;

DROP TABLE rule_conditions;

DROP TABLE rule_goal;

