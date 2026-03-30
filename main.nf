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
    PREFLIGHT(VALIDATE_MANIFEST.out.validated_manifest)
    PROBE_CRAM_REFERENCE(PREFLIGHT.out.context)
    VALIDATE_REFERENCE_BUNDLE(PREFLIGHT.out.context, PROBE_CRAM_REFERENCE.out.reference_probe)
    RUN_MONOGENIC_STAGE(PREFLIGHT.out.context, VALIDATE_REFERENCE_BUNDLE.out.reference_stage)
    RUN_PRS_STAGE(PREFLIGHT.out.context, VALIDATE_REFERENCE_BUNDLE.out.reference_stage)
    RENDER_REPORT(
        PREFLIGHT.out.context,
        PROBE_CRAM_REFERENCE.out.reference_probe,
        VALIDATE_REFERENCE_BUNDLE.out.reference_stage,
        RUN_MONOGENIC_STAGE.out.monogenic_stage,
        RUN_PRS_STAGE.out.prs_stage
    )
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
    PYTHONPATH="${projectDir}" python3 -m workflow.core.validate_manifest \
      --manifest "$manifest" \
      --schema "${projectDir}/manifests/schema/run_manifest.schema.json" \
      --output validated_manifest.json
    """
}

process PREFLIGHT {
    tag "preflight"

    input:
    path validated_manifest

    output:
    path "sample_context.json", emit: context

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.preflight \
      --validated-manifest "$validated_manifest" \
      --output sample_context.json
    """
}

process PROBE_CRAM_REFERENCE {
    tag "reference-probe"

    input:
    path context_json

    output:
    path "reference_probe.json", emit: reference_probe

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.probe_cram_reference \
      --context "$context_json" \
      --output reference_probe.json
    """
}

process VALIDATE_REFERENCE_BUNDLE {
    tag "reference-bundle"

    input:
    path context_json
    path reference_probe_json

    output:
    path "reference_stage.json", emit: reference_stage

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.validate_reference_bundle \
      --context "$context_json" \
      --reference-probe "$reference_probe_json" \
      --output reference_stage.json
    """
}

process RUN_MONOGENIC_STAGE {
    tag "monogenic"

    input:
    path context_json
    path reference_stage_json

    output:
    path "monogenic_stage.json", emit: monogenic_stage

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.run_monogenic_stage \
      --context "$context_json" \
      --reference-stage "$reference_stage_json" \
      --output monogenic_stage.json
    """
}

process RUN_PRS_STAGE {
    tag "prs"

    input:
    path context_json
    path reference_stage_json

    output:
    path "prs_stage.json", emit: prs_stage

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.run_prs_stage \
      --context "$context_json" \
      --reference-stage "$reference_stage_json" \
      --output prs_stage.json
    """
}

process RENDER_REPORT {
    tag "report"

    input:
    path context_json
    path reference_probe_json
    path reference_stage_json
    path monogenic_stage_json
    path prs_stage_json

    output:
    path "pipeline_summary.json", emit: pipeline_summary

    script:
    """
    PYTHONPATH="${projectDir}" python3 -m workflow.core.render_unified_report \
      --context "$context_json" \
      --reference-probe "$reference_probe_json" \
      --reference-stage "$reference_stage_json" \
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
    PYTHONPATH="${projectDir}" python3 -m workflow.core.enforce_pipeline_outcome \
      --summary "$pipeline_summary_json"
    """
}
