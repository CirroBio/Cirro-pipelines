
# Cirro-pipelines

Pipeline configurations for Cirro


## Files

| File name               | Description                                                                                                    | Documentation |
| ----------------------- | -------------------------------------------------------------------------------------------------------------- | ------------- |
| process-definition.json | Database entry for registering the process                                                                     | [schemas/process-definition.schema.json](./schemas/process-definition.schema.json) |
| process-form.json       | Renders the run analysis form                                                                                  | [Customizing the Form](https://docs.cirro.bio/pipelines/configuring-pipeline/#customizing-the-form) |
| process-input.json      | Maps the workflow parameters to the parameters from Cirro                                                      | [Form-To-Workflow Mapping](https://docs.cirro.bio/pipelines/configuring-pipeline/#form-to-workflow-input-mapping) |
| process-output.json     | Post process tasks, these are typically steps needed to create the web-optimized assets used in visualizations | [Output Files](https://docs.cirro.bio/pipelines/configuring-pipeline/#output-files) |
| process-compute.config  | Nextflow/Cromwell compute environment overrides                                                                | [Workflow Configuration override](https://docs.cirro.bio/pipelines/configuring-pipeline/#workflow-configuration-override) |
| preprocess.py          | Python script to preprocess the input data before running the workflow                                          | [Preprocess](https://docs.cirro.bio/pipelines/preprocess-script/) |
| README.md          | Documentation for the pipeline, displayed in the side panel of the pipeline run form                                          | [Providing Documentation](https://docs.cirro.bio/pipelines/configuring-pipeline/#providing-documentation) |