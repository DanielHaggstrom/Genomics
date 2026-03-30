nextflow.enable.dsl = 2

def resolveManifestPath(String commandLine) {
    if (!commandLine) {
        return null
    }

    def matcher = commandLine =~ /(?:^|\s)-params-file(?:=|\s+)(?:"([^"]+)"|'([^']+)'|(\S+))/
    if (!matcher.find()) {
        return null
    }

    return matcher.group(1) ?: matcher.group(2) ?: matcher.group(3)
}

def manifestPath = params.run_manifest ?: resolveManifestPath(workflow.commandLine)
if (!manifestPath) {
    error "Run the pipeline with -params-file manifests/<sample>.yaml"
}

def manifestFile = file(manifestPath, checkIfExists: true)

workflow {
    Channel.of(manifestFile).set { manifest_ch }

    VALIDATE_MANIFEST(manifest_ch)
    PREPARE_CONTEXT(VALIDATE_MANIFEST.out.validated_manifest)
    RUN_MONOGENIC_STAGE(PREPARE_CONTEXT.out.context)
    RUN_PRS_STAGE(PREPARE_CONTEXT.out.context)
    RENDER_REPORT(PREPARE_CONTEXT.out.context, RUN_MONOGENIC_STAGE.out.monogenic_stage, RUN_PRS_STAGE.out.prs_stage)
    ENFORCE_OUTCOME(RENDER_REPORT.out.pipeline_summary)
}

process VALIDATE_MANIFEST {
    tag "validate-manifest"

    input:
    path manifest

    output:
    path "validated_manifest.json", emit: validated_manifest

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.manifest validate \
      --manifest "$manifest" \
      --schema "${projectDir}/manifests/schema/run_manifest.schema.json" \
      --output validated_manifest.json
    """
}

process PREPARE_CONTEXT {
    tag "prepare-context"

    input:
    path validated_manifest

    output:
    path "prepared_context.json", emit: context

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.reference prepare-context \
      --validated-manifest "$validated_manifest" \
      --output prepared_context.json
    """
}

process RUN_MONOGENIC_STAGE {
    tag "monogenic"

    input:
    path context_json

    output:
    path "monogenic_stage.json", emit: monogenic_stage

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.stages.monogenic \
      --context "$context_json" \
      --output monogenic_stage.json
    """
}

process RUN_PRS_STAGE {
    tag "prs"

    input:
    path context_json

    output:
    path "prs_stage.json", emit: prs_stage

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.stages.prs \
      --context "$context_json" \
      --output prs_stage.json
    """
}

process RENDER_REPORT {
    tag "report"

    input:
    path context_json
    path monogenic_stage_json
    path prs_stage_json

    output:
    path "pipeline_summary.json", emit: pipeline_summary

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.reporting.report render \
      --context "$context_json" \
      --monogenic-stage "$monogenic_stage_json" \
      --prs-stage "$prs_stage_json" \
      --output pipeline_summary.json
    """
}

process ENFORCE_OUTCOME {
    tag "enforce-outcome"

    input:
    path pipeline_summary_json

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.reporting.report enforce \
      --summary "$pipeline_summary_json"
    """
}
