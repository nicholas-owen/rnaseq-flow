#!/usr/bin/env nextflow

/*
 * Unit test for rmatsLibType() in workflows/rnaseq.nf.
 *
 *     nextflow run ./tests/unit_rmats_libtype.nf
 *
 * The ./ is required. A script path containing a slash is read as a
 * <owner>/<repo> project name, so `nextflow run tests/unit_rmats_libtype.nf`
 * goes looking for nextflow-io/tests on GitHub.
 *
 * The strandedness mapping is the one place in this pipeline where getting a
 * value the wrong way round is silent. rMATS accepts any --libType and
 * produces a complete set of PSI values either way, so fr-firststrand where
 * fr-secondstrand belongs is a wrong answer and not a crash. That is what B1
 * was.
 *
 * The cohort-disagreement branch has no cohort to exercise it: all six yeast
 * samples are reverse-stranded and GSE151251 is unstranded. Faking a
 * strandedness file to produce a disagreement would test the fake, but the
 * collapse function itself is pure and takes a plain list, so it can be called
 * directly with no channels, no containers and no data.
 *
 * Every call is made in the workflow body rather than in a helper: a function
 * reached through an include resolves there, and a local def calling it is the
 * closure-as-function trap this codebase has walked into before.
 *
 * The run prints six "no process matching config selector" warnings, for the
 * index builders and the two download processes. Those selectors are in
 * conf/modules.config and this script includes no process at all, so the
 * warnings are expected.
 */

include { rmatsLibType } from '../workflows/rnaseq.nf'

workflow {
    def failures = []

    // The dUTP convention. reverse -> fr-firststrand is the assertion that
    // matters: it is the yeast example's own library type, and the transposed
    // value would have run to completion.
    def got_unstranded = rmatsLibType([ [ 'A', 'unstranded' ], [ 'B', 'unstranded' ] ])
    if (got_unstranded != 'fr-unstranded') {
        failures << "unstranded gave '${got_unstranded}', wanted 'fr-unstranded'"
    }

    def got_forward = rmatsLibType([ [ 'A', 'forward' ], [ 'B', 'forward' ] ])
    if (got_forward != 'fr-secondstrand') {
        failures << "forward gave '${got_forward}', wanted 'fr-secondstrand'"
    }

    def got_reverse = rmatsLibType([ [ 'A', 'reverse' ], [ 'B', 'reverse' ] ])
    if (got_reverse != 'fr-firststrand') {
        failures << "reverse gave '${got_reverse}', wanted 'fr-firststrand'"
    }

    // RSEQC writes its verdict with echo, and the workflow reads the file with
    // .text, so the value arrives with a trailing newline.
    def got_padded = rmatsLibType([ [ 'A', 'reverse\n' ], [ 'B', ' reverse ' ] ])
    if (got_padded != 'fr-firststrand') {
        failures << "a verdict with surrounding whitespace gave '${got_padded}', wanted 'fr-firststrand'"
    }

    // A cohort that disagrees has to stop, naming every sample in sample order,
    // which is what lets the offending one be found without opening six files.
    //
    // The message is read from the cause rather than from the exception.
    // error() throws WorkflowScriptErrorException(message), but a function
    // reached through an include is invoked reflectively, so what arrives here
    // is a java.lang.reflect.InvocationTargetException whose own message is
    // null and whose cause carries the text. The ?: covers both, and the class
    // is logged so a change in that wrapping shows up rather than passing
    // silently.
    def mixed_stopped = false
    def mixed_message = null
    try {
        rmatsLibType([ [ 'B', 'reverse' ], [ 'A', 'forward' ], [ 'C', 'reverse' ] ])
    }
    catch (Exception e) {
        mixed_stopped = true
        mixed_message = e.message ?: e.cause?.message
        log.info "  a disagreeing cohort raised ${e.getClass().getName()}: ${mixed_message}"
    }
    if (!mixed_stopped) {
        failures << "a cohort of forward and reverse samples was accepted"
    }
    else if (!mixed_message?.contains('A=forward, B=reverse, C=reverse')) {
        failures << "a disagreeing cohort did not name every sample in order: ${mixed_message}"
    }

    // An unmappable verdict stops rather than returning null, which would
    // reach the command line as --libType null.
    def unknown_stopped = false
    try {
        rmatsLibType([ [ 'A', 'ISR' ] ])
    }
    catch (Exception ignored) {
        unknown_stopped = true
    }
    if (!unknown_stopped) {
        failures << "the Salmon library code 'ISR' was accepted as an RSeQC verdict"
    }

    // An empty verdict is the shape a touch-only RSEQC stub would produce.
    def empty_stopped = false
    try {
        rmatsLibType([ [ 'A', '' ] ])
    }
    catch (Exception ignored) {
        empty_stopped = true
    }
    if (!empty_stopped) {
        failures << "an empty strandedness verdict was accepted"
    }

    // Logged before the error, not inside it: the message handed to error()
    // does not reliably survive, which is the same thing the disagreement
    // assertion above works around.
    if (failures) {
        failures.each { f -> log.error "  FAIL  ${f}" }
        error "rmatsLibType() unit test failed, ${failures.size()} of 7 assertion(s). See the FAIL lines above."
    }
    log.info "rmatsLibType(): 7 assertions passed."
}
