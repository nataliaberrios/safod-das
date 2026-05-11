from pathlib import Path
from textwrap import wrap
from PIL import Image, ImageDraw, ImageFont


ROOT = Path("/Users/nberrios/Documents/Claude/Projects/DAS")
REPO = ROOT / "safod-das"
RESULTS = ROOT / "sanity_results"
OUT = REPO / "reports" / "SAFOD_DAS_session_handoff_2026-05-11.pdf"

PAGE_W, PAGE_H = 2550, 3300
MARGIN_X = 190
MARGIN_TOP = 180
MARGIN_BOTTOM = 170
ACCENT = (28, 80, 128)
MUTED = (92, 100, 112)
TEXT = (30, 36, 42)
LIGHT = (241, 245, 249)
LINE = (198, 208, 220)

FONT = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_ITALIC = "/System/Library/Fonts/Supplemental/Arial Italic.ttf"
FONT_MONO = "/System/Library/Fonts/SFNSMono.ttf"


def font(path, size):
    return ImageFont.truetype(path, size=size)


F_TITLE = font(FONT_BOLD, 70)
F_SUBTITLE = font(FONT, 34)
F_H1 = font(FONT_BOLD, 42)
F_H2 = font(FONT_BOLD, 32)
F_BODY = font(FONT, 27)
F_BODY_BOLD = font(FONT_BOLD, 27)
F_SMALL = font(FONT, 22)
F_SMALL_BOLD = font(FONT_BOLD, 22)
F_MONO = font(FONT_MONO, 21)


def text_height(draw, text, fnt):
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[3] - box[1]


def fit_text(draw, text, fnt, width):
    words = text.split()
    lines = []
    cur = ""
    for word in words:
        test = word if not cur else f"{cur} {word}"
        if draw.textlength(test, font=fnt) <= width:
            cur = test
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines or [""]


class Report:
    def __init__(self):
        self.pages = []
        self.page = None
        self.draw = None
        self.y = MARGIN_TOP
        self.page_no = 0
        self.new_page()

    def new_page(self):
        if self.page is not None:
            self.footer()
            self.pages.append(self.page)
        self.page_no += 1
        self.page = Image.new("RGB", (PAGE_W, PAGE_H), "white")
        self.draw = ImageDraw.Draw(self.page)
        self.y = MARGIN_TOP

    def footer(self):
        y = PAGE_H - 105
        self.draw.line((MARGIN_X, y, PAGE_W - MARGIN_X, y), fill=LINE, width=2)
        self.draw.text((MARGIN_X, y + 26), "SAFOD DAS sanity/reproducibility handoff", font=F_SMALL, fill=MUTED)
        page_text = f"Page {self.page_no}"
        self.draw.text((PAGE_W - MARGIN_X - self.draw.textlength(page_text, font=F_SMALL), y + 26), page_text, font=F_SMALL, fill=MUTED)

    def ensure(self, needed):
        if self.y + needed > PAGE_H - MARGIN_BOTTOM:
            self.new_page()

    def title_page(self):
        self.draw.rectangle((0, 0, PAGE_W, 420), fill=(232, 240, 248))
        self.draw.rectangle((0, 0, 26, PAGE_H), fill=ACCENT)
        self.draw.text((MARGIN_X, 150), "SAFOD DAS Session Handoff", font=F_TITLE, fill=TEXT)
        self.draw.text((MARGIN_X, 250), "What changed after the preflight output, what we learned, and what to run next", font=F_SUBTITLE, fill=MUTED)
        self.draw.text((MARGIN_X, 355), "Generated: 2026-05-11", font=F_SMALL_BOLD, fill=ACCENT)
        self.y = 540
        self.callout(
            "Executive read",
            "The first one-day run was a useful diagnostic, but not a final verdict on the 2024 data. It used daytime-only files, which is not stated in Lellouch et al. 2019. The next controlled tests are all-hours stacking, multi-day stacking, source-channel sweeps, and a same-code test on the 2017 data.",
        )
        self.section("Current State In One Sentence")
        self.paragraph(
            "We have a cleaner SAFOD sanity pipeline, the one-day daytime CC did not reproduce a clear Fig. 7c-style body-wave moveout, and the immediate next work is to test whether that changes when using all continuous files, more days, different virtual-source channels, and finally the 2017 data."
        )
        self.section("Important Correction")
        self.bullets([
            "Gauge length is not the likely blocker in the 5-20 Hz band. The finite-gauge sinc term is close to flat for both 3200 m/s body waves and ~1500 m/s guided/tube waves.",
            "The daytime-only file selection was our hypothesis about cultural noise, not something Lellouch explicitly described.",
            "The vertical fiber geometry makes standard lateral surface-wave hydrology dv/v a poor fit, but does not rule out body-wave, coda, repeating-earthquake, or carefully interpreted guided-wave monitoring.",
        ])

    def section(self, title):
        self.ensure(95)
        if self.y > MARGIN_TOP + 20:
            self.y += 22
        self.draw.text((MARGIN_X, self.y), title, font=F_H1, fill=ACCENT)
        self.y += 58
        self.draw.line((MARGIN_X, self.y, PAGE_W - MARGIN_X, self.y), fill=LINE, width=2)
        self.y += 28

    def subsection(self, title):
        self.ensure(70)
        self.y += 12
        self.draw.text((MARGIN_X, self.y), title, font=F_H2, fill=TEXT)
        self.y += 48

    def paragraph(self, text, fnt=F_BODY, color=TEXT, spacing=16):
        lines = fit_text(self.draw, text, fnt, PAGE_W - 2 * MARGIN_X)
        self.ensure(len(lines) * 36 + spacing)
        for line in lines:
            self.draw.text((MARGIN_X, self.y), line, font=fnt, fill=color)
            self.y += 36
        self.y += spacing

    def bullets(self, items):
        for item in items:
            lines = fit_text(self.draw, item, F_BODY, PAGE_W - 2 * MARGIN_X - 70)
            self.ensure(len(lines) * 36 + 22)
            self.draw.ellipse((MARGIN_X + 8, self.y + 11, MARGIN_X + 22, self.y + 25), fill=ACCENT)
            for i, line in enumerate(lines):
                self.draw.text((MARGIN_X + 55, self.y), line, font=F_BODY, fill=TEXT)
                self.y += 36
            self.y += 12
        self.y += 6

    def numbered(self, items):
        for idx, item in enumerate(items, 1):
            lines = fit_text(self.draw, item, F_BODY, PAGE_W - 2 * MARGIN_X - 90)
            self.ensure(len(lines) * 36 + 22)
            label = f"{idx}."
            self.draw.text((MARGIN_X, self.y), label, font=F_BODY_BOLD, fill=ACCENT)
            for line in lines:
                self.draw.text((MARGIN_X + 70, self.y), line, font=F_BODY, fill=TEXT)
                self.y += 36
            self.y += 12
        self.y += 6

    def code(self, text):
        lines = []
        for raw in text.strip("\n").splitlines():
            if not raw:
                lines.append("")
                continue
            lines.extend(wrap(raw, width=106, break_long_words=False, break_on_hyphens=False) or [""])
        needed = len(lines) * 30 + 42
        self.ensure(needed)
        x0, y0 = MARGIN_X, self.y
        x1 = PAGE_W - MARGIN_X
        y1 = self.y + needed
        self.draw.rounded_rectangle((x0, y0, x1, y1), radius=16, fill=(247, 249, 252), outline=LINE, width=2)
        self.y += 22
        for line in lines:
            self.draw.text((MARGIN_X + 24, self.y), line, font=F_MONO, fill=(26, 43, 60))
            self.y += 30
        self.y += 30

    def callout(self, label, text):
        lines = fit_text(self.draw, text, F_BODY, PAGE_W - 2 * MARGIN_X - 70)
        needed = 72 + len(lines) * 36 + 30
        self.ensure(needed)
        x0, y0 = MARGIN_X, self.y
        x1, y1 = PAGE_W - MARGIN_X, self.y + needed
        self.draw.rounded_rectangle((x0, y0, x1, y1), radius=18, fill=LIGHT, outline=(195, 210, 226), width=2)
        self.draw.rectangle((x0, y0, x0 + 18, y1), fill=ACCENT)
        self.draw.text((x0 + 45, y0 + 24), label, font=F_SMALL_BOLD, fill=ACCENT)
        yy = y0 + 68
        for line in lines:
            self.draw.text((x0 + 45, yy), line, font=F_BODY, fill=TEXT)
            yy += 36
        self.y = y1 + 26

    def image(self, path, caption, max_h=1120):
        path = Path(path)
        if not path.exists():
            self.paragraph(f"[Missing image: {path}]", color=(150, 50, 50))
            return
        img = Image.open(path).convert("RGB")
        max_w = PAGE_W - 2 * MARGIN_X
        scale = min(max_w / img.width, max_h / img.height)
        new_size = (int(img.width * scale), int(img.height * scale))
        needed = new_size[1] + 90
        self.ensure(needed)
        x = (PAGE_W - new_size[0]) // 2
        self.page.paste(img.resize(new_size, Image.LANCZOS), (x, self.y))
        self.y += new_size[1] + 18
        lines = fit_text(self.draw, caption, F_SMALL, max_w)
        for line in lines:
            self.draw.text((MARGIN_X, self.y), line, font=F_SMALL, fill=MUTED)
            self.y += 29
        self.y += 26

    def table(self, headers, rows, widths):
        table_w = PAGE_W - 2 * MARGIN_X
        col_w = [int(table_w * w) for w in widths]
        col_w[-1] = table_w - sum(col_w[:-1])
        row_pad = 16
        x0 = MARGIN_X

        def row_height(cells, fonts):
            max_lines = 1
            for cell, width, fnt in zip(cells, col_w, fonts):
                max_lines = max(max_lines, len(fit_text(self.draw, cell, fnt, width - 28)))
            return max(64, max_lines * 31 + row_pad * 2)

        for row_idx, cells in enumerate([headers] + rows):
            fonts = [F_SMALL_BOLD] * len(headers) if row_idx == 0 else [F_SMALL] * len(headers)
            h = row_height(cells, fonts)
            self.ensure(h + 4)
            fill = (232, 240, 248) if row_idx == 0 else (255, 255, 255)
            x = x0
            for col_idx, (cell, width, fnt) in enumerate(zip(cells, col_w, fonts)):
                self.draw.rectangle((x, self.y, x + width, self.y + h), fill=fill, outline=LINE, width=2)
                yy = self.y + row_pad
                for line in fit_text(self.draw, cell, fnt, width - 28):
                    self.draw.text((x + 14, yy), line, font=fnt, fill=TEXT if row_idx else ACCENT)
                    yy += 31
                x += width
            self.y += h
        self.y += 28

    def finish(self):
        self.footer()
        self.pages.append(self.page)
        OUT.parent.mkdir(parents=True, exist_ok=True)
        self.pages[0].save(OUT, "PDF", resolution=150, save_all=True, append_images=self.pages[1:])


def build():
    r = Report()
    r.title_page()

    r.new_page()
    r.section("Timeline Since Preflight")
    r.table(
        ["Step", "What Happened", "Interpretation"],
        [
            ["Preflight outputs inspected", "Looked at RMS, time-distance, and F-K plots for 2024-10-23 20-21 UTC, channels 150-800.", "Preflight justified continuing to CC, but it was not by itself proof of a recovered body-wave Green's function."],
            ["One-day CC completed", "SLURM job 24589057 ran the CC stage for 2024-10-23 with ch 150-800, source ch 150, 5-20 Hz, 30 s windows, 50% overlap, running-AM 0.1 s, whitening off.", "The run selected 600 daytime continuous files and produced n_stack=1800. Important caveat: daytime-only was our assumption, not stated in Lellouch."],
            ["One-day plots reviewed", "Raw and R+/-10 pre-shifted stacked plots were dominated by near-zero-lag vertical energy, with no clean arrival tracking the 3200 m/s line.", "This did not reproduce Lellouch Fig. 7c, but it is only a daytime one-day result, not a final statement about 2024."],
            ["Follow-up scope corrected", "Decided not to pivot after one failed run. Added tests for more days, alternate virtual-source channels, and same-code testing on 2017 data.", "This is the right scientific escalation: rule out implementation and selection effects before interpreting 2024 as different."],
            ["Hour selection corrected", "Recognized that Lellouch described one day of continuous data, not daytime-only. Added a local code change for SANITY_HOURS=all|day|night, default all.", "The strict reproduction path should use all continuous files. Day/night become explicit source-distribution tests."],
        ],
        [0.18, 0.46, 0.36],
    )

    r.section("Preflight Evidence")
    r.paragraph("The preflight diagnostics were useful for deciding whether the run was worth attempting, but they should not be over-read as a successful Green's-function recovery. They characterize the raw wavefield, not the final stacked correlogram.")
    r.image(RESULTS / "preflight_out" / "rms_per_channel.png", "Preflight RMS per channel. The one-hour diagnostic showed a prominent high-RMS channel near ~310. The later CC job's automated first-files RMS pass flagged channels 705, 715, 721, and 731, so bad-channel behavior appears file/time dependent.", max_h=760)
    r.image(RESULTS / "preflight_out" / "fk_one_file.png", "Preflight F-K spectrum. Energy exists in the 5-20 Hz band and is not purely surface-wave dominated, but the spectrum is broad enough that CC/post-processing is still the decisive test.", max_h=980)
    r.image(RESULTS / "preflight_out" / "timedist_one_file.png", "Raw time-distance panel for the same preflight hour. Useful for spotting coherent bands, bad traces, and common-mode structure before committing cluster time.", max_h=850)

    r.new_page()
    r.section("One-Day CC Result")
    r.paragraph("The one-day CC completed cleanly. The key log facts were: 600 daytime continuous files selected, fs=62.5 Hz, ch 150-800, source ch 150, n_stack=1800, whitening off, and output saved to the sanity_v1 directory. The job warning from temporal normalization was handled with nan_to_num and was not a crash.")
    r.image(RESULTS / "plot_out" / "cc_stacked_2024-10-23.png", "R+/-10 pre-shifted stacked CC for 2024-10-23. Dominant near-zero-lag vertical energy is visible; no clean coherent arrival follows the +/-3200 m/s reference moveout.", max_h=1050)
    r.image(RESULTS / "plot_out" / "cc_fk_2024-10-23.png", "F-K of stacked CC. The result is broad/smeared, with strong near-zero-wavenumber/near-zero-lag structure rather than a clean 3200 m/s ridge.", max_h=930)
    r.callout("Conclusion from this run", "The one-day daytime-only run did not reproduce Lellouch Fig. 7c. The careful conclusion is not 'the 2024 data cannot work'; it is 'this particular daytime one-day source-channel-150 run did not recover the target body-wave moveout.'")

    r.new_page()
    r.section("Scientific Clarifications")
    r.subsection("Gauge Length")
    r.paragraph("The corrected gauge-length model is H(k_z) = i k_z * sinc(k_z L / 2), where k_z = 2*pi*f/v_app. In the 5-20 Hz band, with L about 16 m, the sinc term is near 0.98 for a 20 Hz, 3200 m/s body-wave-like arrival and near 0.93 for a 20 Hz, 1500 m/s guided/tube-like arrival. That is a small effect. It slightly suppresses slower modes more; it does not amplify tube waves.")
    r.subsection("F-K Filtering")
    r.paragraph("F-K analysis maps energy by temporal frequency and spatial wavenumber along the vertical fiber. It can diagnose or isolate apparent velocities already present in the data. It cannot create a missing body-wave Green's function. It may later be useful for separating modes, but it should not be used to quietly redefine the strict reproduction test.")
    r.subsection("Surface-Wave Hydrology dv/v")
    r.paragraph("The advisor's geometry concern is mostly right for standard ambient-noise surface-wave hydrology dv/v: a vertical borehole fiber does not sample lateral Rayleigh/Love paths the way a surface array or horizontal DAS line does. Hydrology-related monitoring could still be possible with body waves, local-event coda, repeating earthquakes, or guided/tube waves, but the sensitivity and interpretation would be different.")

    r.section("Code And Git Changes")
    r.table(
        ["Area", "Change", "Status"],
        [
            ["Multi-day/source-channel support", "Added SANITY_DATES and SANITY_SOURCE_CH to the sanity CC pipeline; plot filenames now include channel/source tags to avoid collisions.", "Pulled on Sherlock via GitHub update 7d8e9d9."],
            ["Gauge-length docs", "Updated README language to stop claiming 16 m gauge length amplifies tube/casing modes in the 5-20 Hz CC band.", "Committed/pulled in the same update."],
            ["All-hours support", "Locally added SANITY_HOURS=all|day|night with default all; output filenames include hoursall/hoursday/hoursnight.", "Local laptop change at report time; needs commit/push/pull before Sherlock can use it."],
            ["Sherlock runtime path", "Sherlock needs /home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python because the new clone lacks DAS-utilities/python.", "Manually restored on Sherlock from /tmp/run_sanity.sbatch.before_pull; should be committed intentionally later."],
        ],
        [0.25, 0.52, 0.23],
    )

    r.new_page()
    r.section("Sherlock State")
    r.bullets([
        "The new Sherlock clone is /home/groups/ettore88/nberrios/safod-das.",
        "The legacy repository /home/groups/ettore88/nberrios/safod_das_git still supplies DAS-utilities/python.",
        "A stash named 'preserve sherlock run_sanity edit' exists from protecting the local run_sanity.sbatch edit before pulling.",
        "Untracked sanity/plot_out/ and sanity/preflight_out/ directories exist on Sherlock; they are outputs and should not be added to Git.",
        "The one-day daytime CC output and plots are local evidence, but the all-hours code was not yet active on Sherlock at the time the week-long job was submitted.",
    ])
    r.subsection("Running Job")
    r.paragraph("Job 24597004 was submitted as a week-long stack using SANITY_DATES=2024-10-22:2024-10-28, ch 150-800, source ch 150. Because the all-hours code had not yet been pulled to Sherlock, this job should be labeled as the week-long daytime-only stack.")
    r.code("""
# Check status/progress
cd /home/groups/ettore88/nberrios/safod-das/sanity
squeue -u nberrios
tail -20 logs/safod_sanity_24597004.err

# When CC finishes, plot the most recent npz produced by the current code
sbatch --export=ALL,STAGE=plot run_sanity.sbatch
""")

    r.section("Recommended Next Tests")
    r.numbered([
        "Let job 24597004 finish. Treat it as the week-long daytime-only cultural-noise stack, not the Lellouch-faithful all-hours run.",
        "Commit and push the local SANITY_HOURS change from the laptop, then pull it on Sherlock. Confirm the Sherlock run_sanity.sbatch still includes the legacy DAS-utilities path.",
        "Run the all-hours week stack with SANITY_HOURS=all. This is the fairer Lellouch-style multi-day test.",
        "Run day and night week stacks only as source-timing comparisons after the all-hours result exists.",
        "Run source-channel sweeps for at least source channels 150, 180, 200, and 250 under the same all-hours setup.",
        "Run the same code on the 2017 Lellouch-era data once its path/manifest is known. If 2017 reproduces and 2024 does not, the result becomes a strong 2017-vs-2024 source/fiber/wavefield comparison.",
    ])

    r.new_page()
    r.section("Exact Commands For The Next Clean All-Hours Run")
    r.paragraph("Use these only after the SANITY_HOURS code has been committed, pushed, pulled on Sherlock, and the legacy DAS-utilities path has been preserved in the active run_sanity.sbatch.")
    r.code("""
cd /home/groups/ettore88/nberrios/safod-das/sanity

export SANITY_DATES=2024-10-22:2024-10-28
export SANITY_HOURS=all
export SANITY_CH_START=150
export SANITY_CH_END=800
export SANITY_SOURCE_CH=150

sbatch --export=ALL,STAGE=cc run_sanity.sbatch
""")
    r.subsection("Expected Output Naming")
    r.code("""
/oak/stanford/groups/ettore88/nberrios/sanity_v1/
  sanity_cc_2024-10-22_to_2024-10-28_7d_hoursall_ch150-800_src150.npz

/home/groups/ettore88/nberrios/safod-das/sanity/plot_out/
  cc_raw_2024-10-22_to_2024-10-28_7d_hoursall_ch150-800_src150.png
  cc_stacked_2024-10-22_to_2024-10-28_7d_hoursall_ch150-800_src150.png
  cc_fk_2024-10-22_to_2024-10-28_7d_hoursall_ch150-800_src150.png
""")
    r.callout("Decision rule", "Do not treat one failed daytime run as the final outcome. The stronger decision point is: all-hours multi-day + source-channel sweep + 2017 same-code test.")

    r.finish()
    print(OUT)


if __name__ == "__main__":
    build()
