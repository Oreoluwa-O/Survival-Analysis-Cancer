import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

# ── 1. Load data ───
# Using a simulated TCGA-style dataset for reproducibility.
# Replace this section with actual TCGA clinical data download from:
# https://portal.gdc.cancer.gov/

np.random.seed(42)
n = 300

data = pd.DataFrame({
    "duration":    np.random.exponential(scale=24, size=n).round(1),  # months
    "event":       np.random.binomial(1, 0.65, size=n),               # 1 = death observed
    "stage":       np.random.choice(["Stage I", "Stage II",
                                     "Stage III", "Stage IV"], size=n,
                                    p=[0.25, 0.30, 0.25, 0.20]),
    "age":         np.random.randint(40, 80, size=n),
    "gender":      np.random.choice(["Male", "Female"], size=n),
})

print("Dataset shape:", data.shape)
print(data.head())

# ── 2. Kaplan-Meier survival curves by cancer stage ───
fig, ax = plt.subplots(figsize=(10, 6))

for stage in ["Stage I", "Stage II", "Stage III", "Stage IV"]:
    mask = data["stage"] == stage
    kmf = KaplanMeierFitter()
    kmf.fit(
        durations=data.loc[mask, "duration"],
        event_observed=data.loc[mask, "event"],
        label=stage
    )
    kmf.plot_survival_function(ax=ax, ci_show=True)

ax.set_title("Kaplan-Meier Survival Curves by Cancer Stage", fontsize=14)
ax.set_xlabel("Time (months)")
ax.set_ylabel("Survival probability")
ax.legend(title="Cancer Stage")
plt.tight_layout()
plt.savefig("outputs/kaplan_meier_by_stage.png", dpi=150)
plt.show()
print("Saved: outputs/kaplan_meier_by_stage.png")

# ── 3. Log-rank test: Stage I vs Stage IV ──
s1 = data[data["stage"] == "Stage I"]
s4 = data[data["stage"] == "Stage IV"]

result = logrank_test(
    s1["duration"], s4["duration"],
    event_observed_A=s1["event"],
    event_observed_B=s4["event"]
)

print("\nLog-rank test — Stage I vs Stage IV:")
print(f"  Test statistic : {result.test_statistic:.4f}")
print(f"  p-value        : {result.p_value:.4f}")
if result.p_value < 0.05:
    print("  Conclusion: Survival curves differ significantly (p < 0.05)")
else:
    print("  Conclusion: No significant difference detected")

# ── 4. Cox Proportional Hazards model ───
Cox_data = data.copy()
Cox_data["stage_num"] = Cox_data["stage"].map(
    {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4}
)
Cox_data["gender_num"] = (Cox_data["gender"] == "Male").astype(int)

CPH = CoxPHFitter()
CPH.fit(
    Cox_data[["duration", "event", "stage_num", "age", "gender_num"]],
    duration_col="duration",
    event_col="event"
)

print("\nCox Proportional Hazards Model Summary:")
CPH.print_summary()

fig2, ax2 = plt.subplots(figsize=(8, 4))
CPH.plot(ax=ax2)
ax2.set_title("Cox Model — Hazard Ratios with 95% Confidence Intervals")
plt.tight_layout()
plt.savefig("outputs/cox_hazard_ratios.png", dpi=150)
plt.show()
print("Saved: outputs/cox_hazard_ratios.png")
