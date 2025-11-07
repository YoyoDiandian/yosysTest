#!/bin/bash

# Get script directory
DIRNAME=$0
if [ "${DIRNAME:0:1}" = "/" ];then
    CURDIR=`dirname $DIRNAME`
else
    CURDIR="`pwd`"/"`dirname $DIRNAME`"
fi

echo "Starting batch testing..."
echo "Working directory: $CURDIR"

# Demo directory containing input verilog files
DEMO_DIR="$CURDIR/techlibs/pango/demo"
OUTPUT_DIR="$CURDIR/techlibs/pango/outputs"
ERROR_LOG="$OUTPUT_DIR/batch_test_errors.log"

# Timeout in seconds (30 minutes = 1800 seconds)
TIMEOUT=1800

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Clear previous error log
> "$ERROR_LOG"

# Counter for statistics
TOTAL_FILES=0
SUCCESS_FILES=0
FAILED_FILES=0
TIMEOUT_FILES=0

echo "========================================" | tee -a "$ERROR_LOG"
echo "Batch Test Started: $(date)" | tee -a "$ERROR_LOG"
echo "Timeout: ${TIMEOUT} seconds (30 minutes)" | tee -a "$ERROR_LOG"
echo "========================================" | tee -a "$ERROR_LOG"
echo ""

# Find all .v files in demo directory
for vfile in "$DEMO_DIR"/*.v; do
    # Skip if no .v files found
    if [ ! -f "$vfile" ]; then
        echo "No .v files found in $DEMO_DIR"
        exit 1
    fi
    
    TOTAL_FILES=$((TOTAL_FILES + 1))
    filename=$(basename "$vfile")
    
    echo "----------------------------------------"
    echo "[$TOTAL_FILES] Processing: $filename"
    echo "Started at: $(date)"
    echo "----------------------------------------"
    
    # Record start time
    START_TIME=$(date +%s)
    
    # Run test.sh with timeout (macOS compatible)
    # Run the test in background and monitor it
    bash "$CURDIR/test.sh" "$vfile" > "$OUTPUT_DIR/${filename%.v}_test.log" 2>&1 &
    TEST_PID=$!
    
    # Wait for the process with timeout
    WAIT_TIME=0
    TIMED_OUT=0
    while kill -0 $TEST_PID 2>/dev/null; do
        sleep 1
        WAIT_TIME=$((WAIT_TIME + 1))
        if [ $WAIT_TIME -ge $TIMEOUT ]; then
            # Timeout reached, kill the process
            kill -9 $TEST_PID 2>/dev/null
            TIMED_OUT=1
            break
        fi
    done
    
    # Get exit code if process finished normally
    wait $TEST_PID 2>/dev/null
    EXIT_CODE=$?
    
    # Record end time
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    ELAPSED_MIN=$((ELAPSED / 60))
    ELAPSED_SEC=$((ELAPSED % 60))
    
    # Check exit status
    if [ $TIMED_OUT -eq 1 ]; then
        # Timeout occurred
        TIMEOUT_FILES=$((TIMEOUT_FILES + 1))
        FAILED_FILES=$((FAILED_FILES + 1))
        echo "✗ TIMEOUT: $filename exceeded ${TIMEOUT}s limit" | tee -a "$ERROR_LOG"
        echo "  Terminated at: $(date)" | tee -a "$ERROR_LOG"
        echo "  Log file: $OUTPUT_DIR/${filename%.v}_test.log" | tee -a "$ERROR_LOG"
        echo "" >> "$ERROR_LOG"
    elif [ $EXIT_CODE -eq 0 ]; then
        SUCCESS_FILES=$((SUCCESS_FILES + 1))
        echo "✓ SUCCESS: $filename (${ELAPSED_MIN}m ${ELAPSED_SEC}s)"
        echo "[$(date)] SUCCESS: $filename (${ELAPSED_MIN}m ${ELAPSED_SEC}s)" >> "$ERROR_LOG"
    else
        FAILED_FILES=$((FAILED_FILES + 1))
        echo "✗ FAILED: $filename (exit code: $EXIT_CODE, ${ELAPSED_MIN}m ${ELAPSED_SEC}s)" | tee -a "$ERROR_LOG"
        echo "  Log file: $OUTPUT_DIR/${filename%.v}_test.log" | tee -a "$ERROR_LOG"
        echo "" >> "$ERROR_LOG"
    fi
    
    echo ""
done

# Print summary
echo "========================================"
echo "Batch Test Completed: $(date)"
echo "========================================"
echo "Total files:    $TOTAL_FILES"
echo "Successful:     $SUCCESS_FILES"
echo "Failed:         $FAILED_FILES"
echo "  - Timeout:    $TIMEOUT_FILES"
echo "========================================"
echo ""
echo "Error log: $ERROR_LOG"
echo ""

# Write summary to error log
echo "" >> "$ERROR_LOG"
echo "========================================" >> "$ERROR_LOG"
echo "Summary" >> "$ERROR_LOG"
echo "========================================" >> "$ERROR_LOG"
echo "Total files:    $TOTAL_FILES" >> "$ERROR_LOG"
echo "Successful:     $SUCCESS_FILES" >> "$ERROR_LOG"
echo "Failed:         $FAILED_FILES" >> "$ERROR_LOG"
echo "  - Timeout:    $TIMEOUT_FILES" >> "$ERROR_LOG"
echo "========================================" >> "$ERROR_LOG"

# Exit with error code if any test failed
if [ $FAILED_FILES -gt 0 ]; then
    exit 1
else
    exit 0
fi
