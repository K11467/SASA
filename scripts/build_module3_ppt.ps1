param(
    [string]$InputPptm = "D:\appData\weiXin\xwechat_files\wxid_0ba8c32vmkmv22_0828\msg\file\2026-06\SASA_PPI_module2_method_pages.pptm",
    [string]$OutputPptm = "D:\courseProject\code\SASAGpu\paper\SASA_PPI_module2_module3_visual.pptx"
)

$ErrorActionPreference = "Stop"

$ppLayoutBlank = 12
$ppSaveAsOpenXMLPresentationMacroEnabled = 25
$ppSaveAsOpenXMLPresentation = 24
$msoShapeRectangle = 1
$msoShapeRoundedRectangle = 5
$msoShapeChevron = 52
$msoShapeOval = 9
$msoShapeLine = 9
$msoTextOrientationHorizontal = 1

$navy = 0x604F0E
$teal = 0x787846
$lightTeal = 0xD9F0F0
$orange = 0x0A5FE9
$gray = 0xF5F3F1
$midGray = 0xA6A6A6
$dark = 0x33261E
$white = 0xFFFFFF
$red = 0x3248C9
$green = 0x527533

function Set-Fill($shape, [int]$color, [double]$transparency = 0) {
    $shape.Fill.Visible = -1
    $shape.Fill.ForeColor.RGB = $color
    $shape.Fill.Transparency = $transparency
}

function Set-Line($shape, [int]$color, [double]$weight = 1) {
    $shape.Line.Visible = -1
    $shape.Line.ForeColor.RGB = $color
    $shape.Line.Weight = $weight
}

function Add-Text(
    $slide,
    [double]$left,
    [double]$top,
    [double]$width,
    [double]$height,
    [string]$text,
    [double]$fontSize = 18,
    [int]$color = $dark,
    [bool]$bold = $false,
    [string]$font = "Aptos",
    [int]$align = 1
) {
    $shape = $slide.Shapes.AddTextbox($msoTextOrientationHorizontal, $left, $top, $width, $height)
    $shape.TextFrame.MarginLeft = 0
    $shape.TextFrame.MarginRight = 0
    $shape.TextFrame.MarginTop = 0
    $shape.TextFrame.MarginBottom = 0
    $shape.TextFrame.WordWrap = -1
    $range = $shape.TextFrame.TextRange
    $range.Text = $text
    $range.Font.Name = $font
    $range.Font.Size = $fontSize
    $range.Font.Color.RGB = $color
    $range.Font.Bold = $(if ($bold) { -1 } else { 0 })
    $range.ParagraphFormat.Alignment = $align
    return $shape
}

function Add-Rect(
    $slide,
    [double]$left,
    [double]$top,
    [double]$width,
    [double]$height,
    [int]$fill,
    [int]$line = $fill,
    [bool]$rounded = $false
) {
    $shapeType = $(if ($rounded) { $msoShapeRoundedRectangle } else { $msoShapeRectangle })
    $shape = $slide.Shapes.AddShape($shapeType, $left, $top, $width, $height)
    Set-Fill $shape $fill
    Set-Line $shape $line 0.8
    return $shape
}

function Add-Header($slide, [string]$section, [string]$title, [int]$page, [int]$total) {
    $section = $section.Replace([string][char]0x8DEF, [string][char]0x00B7)
    $dot = $slide.Shapes.AddShape($msoShapeOval, 43, 43, 15, 15)
    Set-Fill $dot $orange
    $dot.Line.Visible = 0
    Add-Text $slide 64 39 520 24 $section 13 $orange $true "Aptos" | Out-Null
    Add-Text $slide 43 77 870 54 $title 31 $navy $true "Aptos Display" | Out-Null
    Add-Text $slide 900 515 48 18 "$page / $total" 9 $midGray $false "Aptos" 3 | Out-Null
}

function Add-Card(
    $slide,
    [double]$left,
    [double]$top,
    [double]$width,
    [double]$height,
    [string]$kicker,
    [string]$title,
    [string]$body,
    [int]$accent = $orange
) {
    Add-Rect $slide $left $top $width $height $white $lightTeal | Out-Null
    Add-Rect $slide $left $top 8 $height $accent $accent | Out-Null
    Add-Text $slide ($left + 18) ($top + 16) ($width - 32) 18 $kicker 10 $accent $true "Aptos" | Out-Null
    Add-Text $slide ($left + 18) ($top + 42) ($width - 32) 30 $title 17 $navy $true "Aptos Display" | Out-Null
    Add-Text $slide ($left + 18) ($top + 84) ($width - 32) ($height - 94) $body 12 $dark $false "Aptos" | Out-Null
}

function Add-Step($slide, [double]$left, [double]$top, [double]$width, [string]$step, [string]$title, [string]$body, [int]$fill) {
    Add-Rect $slide $left $top $width 116 $fill $fill | Out-Null
    Add-Text $slide ($left + 12) ($top + 16) ($width - 24) 18 $step 10 $lightTeal $true "Aptos" 2 | Out-Null
    Add-Text $slide ($left + 12) ($top + 42) ($width - 24) 32 $title 16 $white $true "Aptos Display" 2 | Out-Null
    Add-Text $slide ($left + 12) ($top + 82) ($width - 24) 26 $body 10 $lightTeal $false "Aptos" 2 | Out-Null
}

function Add-MetricBar($slide, [double]$left, [double]$top, [double]$width, [string]$label, [double]$value, [int]$fill, [string]$suffix = "") {
    Add-Text $slide $left $top 235 20 $label 12 $navy $true "Aptos" | Out-Null
    Add-Rect $slide ($left + 240) ($top + 3) $width 13 $gray $gray $true | Out-Null
    Add-Rect $slide ($left + 240) ($top + 3) ($width * $value) 13 $fill $fill $true | Out-Null
    Add-Text $slide ($left + 250 + $width) ($top - 1) 82 20 ("{0:N4}{1}" -f $value, $suffix) 12 $dark $true "Aptos" | Out-Null
}

function Add-TableCell($slide, [double]$left, [double]$top, [double]$width, [double]$height, [string]$text, [int]$fill, [int]$textColor, [bool]$bold = $false, [int]$align = 1) {
    Add-Rect $slide $left $top $width $height $fill $white | Out-Null
    Add-Text $slide ($left + 6) ($top + 7) ($width - 12) ($height - 10) $text 10 $textColor $bold "Aptos" $align | Out-Null
}

function Add-Line($slide, [double]$x1, [double]$y1, [double]$x2, [double]$y2, [int]$color, [double]$weight = 2) {
    $shape = $slide.Shapes.AddLine($x1, $y1, $x2, $y2)
    Set-Line $shape $color $weight
    return $shape
}

function Add-CircleNode($slide, [double]$left, [double]$top, [double]$size, [string]$label, [string]$body, [int]$fill) {
    $circle = $slide.Shapes.AddShape($msoShapeOval, $left, $top, $size, $size)
    Set-Fill $circle $fill
    $circle.Line.Visible = 0
    Add-Text $slide ($left + 8) ($top + 17) ($size - 16) 24 $label 16 $white $true "Aptos Display" 2 | Out-Null
    Add-Text $slide ($left - 32) ($top + $size + 13) ($size + 64) 42 $body 11 $navy $true "Aptos" 2 | Out-Null
}

function Add-StatCard($slide, [double]$left, [double]$top, [double]$width, [double]$height, [string]$label, [string]$value, [string]$note, [int]$fill) {
    Add-Rect $slide $left $top $width $height $fill $fill $true | Out-Null
    Add-Text $slide ($left + 18) ($top + 16) ($width - 36) 22 $label 11 $lightTeal $true "Aptos" 2 | Out-Null
    Add-Text $slide ($left + 18) ($top + 42) ($width - 36) 48 $value 31 $white $true "Aptos Display" 2 | Out-Null
    Add-Text $slide ($left + 18) ($top + 93) ($width - 36) 38 $note 11 $lightTeal $false "Aptos" 2 | Out-Null
}

if (-not (Test-Path -LiteralPath $InputPptm)) {
    throw "Input PPTM not found: $InputPptm"
}

$outputDir = Split-Path -Parent $OutputPptm
New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

$app = New-Object -ComObject PowerPoint.Application
$app.Visible = -1
$pres = $null

try {
    $pres = $app.Presentations.Open($InputPptm, $false, $false, $false)
    $insertAt = $pres.Slides.Count
    $total = $pres.Slides.Count + 8

    foreach ($slide in $pres.Slides) {
        $hasPageCounter = $false
        foreach ($shape in $slide.Shapes) {
            try {
                if ($shape.HasTextFrame -and $shape.TextFrame.HasText) {
                    $text = $shape.TextFrame.TextRange.Text
                    if ($text -match "^\d+\s*/\s*17$") {
                        $shape.TextFrame.TextRange.Text = "$($slide.SlideIndex) / $total"
                        $hasPageCounter = $true
                    }
                }
            } catch {}
        }
        if (-not $hasPageCounter -and $slide.SlideIndex -ge 2 -and $slide.SlideIndex -lt $pres.Slides.Count) {
            Add-Text $slide 900 515 48 18 "$($slide.SlideIndex) / $total" 9 $midGray $false "Aptos" 3 | Out-Null
        }
    }

    # Slide 17: overview
    $slide = $pres.Slides.Add($insertAt, $ppLayoutBlank)
    Add-Header $slide "04.3  ·  Module 3 Overview" "From weak labels to leakage-controlled EGNN" 17 $total
    $fills = @($teal, $teal, $orange, $navy, $navy)
    $steps = @(
        @("STEP 1", "Weak labels", "Delta-SASA threshold"),
        @("STEP 2", "ESM-2 650M", "1280-D embeddings"),
        @("STEP 3", "Strict features", "No SASA inputs"),
        @("STEP 4", "EGNN", "3D residue graph"),
        @("STEP 5", "Prediction", "Residue probability")
    )
    for ($i = 0; $i -lt 5; $i++) {
        $left = 48 + 182 * $i
        Add-Step $slide $left 175 155 $steps[$i][0] $steps[$i][1] $steps[$i][2] $fills[$i]
        if ($i -lt 4) {
            $chevron = $slide.Shapes.AddShape($msoShapeChevron, ($left + 157), 217, 24, 32)
            Set-Fill $chevron 0xDDD0BE
            $chevron.Line.Visible = 0
        }
    }
    Add-Text $slide 53 355 850 38 "Core design rule: Delta-SASA is the label source, not a shortcut feature." 22 $orange $true "Aptos Display" 2 | Out-Null
    Add-Text $slide 105 411 750 42 "The deployable model learns sequence and geometry while keeping weak-label construction separate from prediction." 15 $dark $false "Aptos" 2 | Out-Null

    # Slide 18: leakage control
    $slide = $pres.Slides.Add($insertAt + 1, $ppLayoutBlank)
    Add-Header $slide "04.3  ·  Leakage-Controlled Features" "Separating label construction from deployable inputs" 18 $total
    Add-Card $slide 52 168 258 242 "LABEL SOURCE" "Delta-SASA" "SASA_apo - SASA_holo`n`nUsed to generate weak interface labels at a 2.0 A^2 threshold." $orange
    Add-Card $slide 351 168 258 242 "STRICT PRIMARY" "ESM + geometry" "1280-D ESM-2 embeddings`nBackbone dihedrals`nHSE + hydrophobicity`n3D coordinates" $green
    Add-Card $slide 650 168 258 242 "DIAGNOSTIC ONLY" "Holo-aware inputs" "Adding SASA_holo strongly raises scores because it is close to the label definition.`n`nDo not report as deployment performance." $red
    Add-Text $slide 82 445 800 36 "Strict model: 1287 scalar node features. Delta-SASA, SASA_holo and SASA_apo are excluded." 15 $navy $true "Aptos" 2 | Out-Null

    # Slide 19: feature representation
    $slide = $pres.Slides.Add($insertAt + 2, $ppLayoutBlank)
    Add-Header $slide "04.4  ·  Residue Representation" "ESM-2 650M semantics meet local structural geometry" 19 $total
    $featureRows = @(
        @("ESM-2 embedding", "1280", "Evolutionary and sequence context", $orange),
        @("Dihedral encoding", "4", "sin(phi), cos(phi), sin(psi), cos(psi)", $teal),
        @("Half-sphere exposure", "2", "Local neighborhood orientation", $teal),
        @("Hydrophobicity", "1", "Residue physicochemical tendency", $teal)
    )
    $y = 162
    foreach ($row in $featureRows) {
        Add-Rect $slide 58 $y 585 64 $white $lightTeal | Out-Null
        Add-Rect $slide 58 $y 8 64 $row[3] $row[3] | Out-Null
        Add-Text $slide 84 ($y + 15) 220 24 $row[0] 15 $navy $true "Aptos Display" | Out-Null
        Add-Text $slide 325 ($y + 15) 72 24 $row[1] 17 $orange $true "Aptos Display" 2 | Out-Null
        Add-Text $slide 420 ($y + 16) 200 24 $row[2] 11 $dark $false "Aptos" | Out-Null
        $y += 72
    }
    Add-Rect $slide 688 177 210 225 $navy $navy $true | Out-Null
    Add-Text $slide 720 206 150 28 "STRICT NODE INPUT" 12 $lightTeal $true "Aptos" 2 | Out-Null
    Add-Text $slide 719 249 150 52 "1287" 42 $white $true "Aptos Display" 2 | Out-Null
    Add-Text $slide 715 309 160 44 "scalar features`nper residue" 15 $lightTeal $false "Aptos" 2 | Out-Null
    Add-Text $slide 695 435 195 36 "Coordinates are used by the graph, not counted in the scalar vector." 11 $dark $false "Aptos" 2 | Out-Null

    # Slide 20: EGNN
    $slide = $pres.Slides.Add($insertAt + 3, $ppLayoutBlank)
    Add-Header $slide "04.4  ·  EGNN Spatial Modeling" "Reasoning over a 3D residue graph" 20 $total
    Add-Card $slide 52 168 258 232 "01  GRAPH" "Residue nodes" "One node per target-chain residue.`nEdges connect C-alpha atoms within an 8 A cutoff." $teal
    Add-Card $slide 351 168 258 232 "02  MESSAGE PASSING" "Equivariant updates" "Messages depend on residue states and normalized pairwise distances.`nCoordinate updates respect spatial symmetry." $orange
    Add-Card $slide 650 168 258 232 "03  OUTPUT" "Interface probability" "Three EGNN layers aggregate local context.`nAn MLP head outputs residue-level interface probabilities." $navy
    Add-Text $slide 121 443 720 30 "Translation, rotation and reflection symmetry are preserved by construction." 17 $navy $true "Aptos Display" 2 | Out-Null

    # Slide 21: leakage ablation
    $slide = $pres.Slides.Add($insertAt + 4, $ppLayoutBlank)
    Add-Header $slide "05  ·  Experiments and Ablation" "Leakage audit changes the interpretation of the result" 21 $total
    Add-MetricBar $slide 73 170 390 "apo-SASA rank" 0.3317 $teal
    Add-MetricBar $slide 73 207 390 "MLP / ESM-only" 0.6290 $teal
    Add-MetricBar $slide 73 244 390 "GCN / ESM + geometry" 0.5961 $teal
    Add-MetricBar $slide 73 281 390 "EGNN / ESM + geometry" 0.7745 $green
    Add-MetricBar $slide 73 318 390 "EGNN / apo + ESM + geometry" 0.7787 $green
    Add-MetricBar $slide 73 355 390 "EGNN / apo+holo + ESM + geometry" 0.8948 $red
    Add-Text $slide 683 139 185 28 "TEST F1" 14 $navy $true "Aptos Display" 2 | Out-Null
    Add-Rect $slide 674 205 202 184 $gray $gray $true | Out-Null
    Add-Text $slide 699 225 152 32 "+0.1203" 34 $red $true "Aptos Display" 2 | Out-Null
    Add-Text $slide 696 273 160 72 "apparent gain after`nadding holo SASA`n`nThis is leakage evidence,`nnot a deployable improvement." 13 $dark $false "Aptos" 2 | Out-Null
    Add-Text $slide 98 442 770 28 "Primary result: strict EGNN reaches F1 = 0.7745, AUROC = 0.9322 and AUPRC = 0.8421." 15 $navy $true "Aptos" 2 | Out-Null

    # Slide 22: external benchmarks
    $slide = $pres.Slides.Add($insertAt + 5, $ppLayoutBlank)
    Add-Header $slide "05  ·  External Benchmarks" "Strict evaluation exposes a real generalization gap" 22 $total
    $columns = @(58, 263, 446, 590, 708, 818)
    $widths = @(205, 183, 144, 118, 110, 90)
    $headers = @("Dataset", "Role", "Model", "F1", "AUROC", "AUPRC")
    for ($i = 0; $i -lt $headers.Count; $i++) {
        Add-TableCell $slide $columns[$i] 165 $widths[$i] 36 $headers[$i] $navy $white $true 2
    }
    $tableRows = @(
        @("Dset_186-local", "Strict primary", "EGNN", "0.3529", "0.7277", "0.3050"),
        @("PDBtest_315-local", "Strict primary", "EGNN", "0.3040", "0.6799", "0.2675"),
        @("PDBtest_315-local", "Diagnostic", "Holo EGNN", "0.6761", "0.8946", "0.7408"),
        @("PDBtest_315-local", "Diagnostic", "Cross-chain", "0.6816", "0.8983", "0.6906")
    )
    for ($rowIndex = 0; $rowIndex -lt $tableRows.Count; $rowIndex++) {
        $fill = $(if ($rowIndex -lt 2) { $white } else { 0xEDE8FA })
        for ($i = 0; $i -lt $headers.Count; $i++) {
            Add-TableCell $slide $columns[$i] (201 + 43 * $rowIndex) $widths[$i] 43 $tableRows[$rowIndex][$i] $fill $dark ($i -eq 1) 2
        }
    }
    Add-Text $slide 76 410 810 54 "Strict external scores are substantially lower than internal scores. The current package is a reproducible baseline, not a state-of-the-art claim." 18 $orange $true "Aptos Display" 2 | Out-Null

    # Slide 23: cross-chain analysis
    $slide = $pres.Slides.Add($insertAt + 6, $ppLayoutBlank)
    Add-Header $slide "05  ·  Cross-chain Analysis" "A useful analysis module, not a stable core improvement" 23 $total
    Add-Card $slide 61 173 255 235 "WHAT IT ADDS" "Partner-chain attention" "Target residues attend to partner-chain C-alpha positions using distance-weighted context." $teal
    Add-Card $slide 352 173 255 235 "WHAT WE OBSERVE" "Metric trade-off" "On holo-aware PDBtest_315-local:`nRecall: 0.6252 -> 0.6629`nF1: 0.6761 -> 0.6816`nAUPRC: 0.7408 -> 0.6906" $orange
    Add-Card $slide 643 173 255 235 "HOW WE CLAIM IT" "Conservative framing" "Cross-chain attention is retained for analysis. Its benefits are benchmark-dependent and not consistently positive." $navy
    Add-Text $slide 99 448 760 26 "Next step: retrain cross-chain variants under strict no-SASA inputs and matched manifests." 15 $navy $true "Aptos" 2 | Out-Null

    # Slide 24: summary
    $slide = $pres.Slides.Add($insertAt + 7, $ppLayoutBlank)
    Add-Header $slide "06  ·  Summary" "What module 3 adds to the complete project" 24 $total
    Add-Card $slide 58 166 270 228 "01  METHOD" "Closed-loop pipeline" "Self-developed SASA`n-> Delta-SASA weak labels`n-> ESM-2 650M + geometry`n-> EGNN prediction" $orange
    Add-Card $slide 345 166 270 228 "02  RIGOR" "Leakage-controlled study" "Strict no-SASA primary model`nApo-only ablation`nHolo-aware diagnostics`nExternal local benchmarks" $green
    Add-Card $slide 632 166 270 228 "03  DELIVERY" "Reproducible package" "Tracked checkpoint`nPrediction CSV exports`nAutomated tests`nWindows GBK compatibility" $navy
    Add-Text $slide 80 438 800 32 "Takeaway: the engineering pipeline is complete; the next research priority is strict external generalization." 16 $orange $true "Aptos Display" 2 | Out-Null

    $pres.SaveAs($OutputPptm, $ppSaveAsOpenXMLPresentationMacroEnabled)
    Write-Output "Saved: $OutputPptm"
    Write-Output "Slides: $($pres.Slides.Count)"
}
finally {
    if ($pres) {
        $pres.Close()
    }
    $app.Quit()
}
