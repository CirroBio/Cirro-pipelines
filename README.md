
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
| preprocess.py           | Python script to preprocess the input data before running the workflow                                         | [Preprocess](https://docs.cirro.bio/pipelines/preprocess-script/) |
| README.md               | Documentation for the pipeline, displayed in the side panel of the pipeline run form                           | [Providing Documentation](https://docs.cirro.bio/pipelines/configuring-pipeline/#providing-documentation) |
## Sizing compute to the input

A `process-compute.config` can choose a task's resources from the size of its input files,
so large samples get large instances while small ones stay cheap. Write the directive as a
closure: Nextflow evaluates it for each task before submitting it, with the process's input
variables in scope.

```groovy
withName: 'BISMARK_ALIGN' {
    cpus   = { (reads*.size().sum() > 10.GB.toBytes() ? 48 : 12) * task.attempt }
    memory = { (reads*.size().sum() > 10.GB.toBytes() ? 192.GB : 72.GB) * task.attempt }
    queue  = { reads*.size().sum() > 10.GB.toBytes() ? env('PW_ONDEMAND_JOB_QUEUE') : env('PW_SPOT_JOB_QUEUE') }
}
```

- `reads` is the name the module gives its input (`tuple val(meta), path(reads)`). Read the
  module's `input:` block for it; a later release that renames the input fails the task
  rather than falling back to the defaults, so the override belongs to one pipeline version.
- `reads*.size().sum()` is the total size in bytes, of one file or of a read pair alike.
  Compare it with `10.GB.toBytes()`, not `10.GB`. The size is of the file as staged
  (usually gzipped), and a step after trimming sees the trimmed reads.
- The smaller tier should repeat the pipeline's own label resources (`process_high` here),
  since this override replaces them. Keep its time limits at least the pipeline's.
- Multiply by `task.attempt` so retries still grow, and set `process.resourceLimits` so a
  late retry is capped at an instance that exists instead of waiting forever.
- Moving the large tier to `PW_ONDEMAND_JOB_QUEUE` keeps multi-day tasks from being lost to
  Spot reclaims.
- Nextflow's strict config syntax allows no helper functions or variables, so the
  expression is repeated in each directive. Use `env('...')`, not `"${...}"`, for queues.

See `nf-core/methylseq/4.2/process-compute.config`.
