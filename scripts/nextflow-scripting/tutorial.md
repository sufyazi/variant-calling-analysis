# Porting a Workflow into Nextflow Pipeline

This will be my attempt at porting the BCFtools workflow into a Nextflow pipeline. I will be using the [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) syntax.

## Step 1: Defining the process
First we will need to define a process to run in the test pipeline. We can first do a simple process that will run the `bcftools` command on a single file. We will call this process `bcftools_process`.

```nextflow