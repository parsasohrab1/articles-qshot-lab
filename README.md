# articles-qshot-lab
## SRS — آزمایشگاه بودجه‌ی شات (Shot-Budget Lab)

### 1. Introduction
- **Purpose:** تعریف یک چارچوب بازتولیدپذیر برای مقایسه‌ی سه استراتژی تخصیص بودجه‌ی اندازه‌گیری (یکنواخت، وزن‌دار، سایه‌های کلاسیک) در تخمین مقدار انتظاری مشاهده‌پذیرهای مجموع پائولی.
- **Scope:** نرم‌افزار شبیه‌سازی برای سیستم‌های کوچک کیوبیتی (۳ تا ۸ کیوبیت)، با خروجی RMSE و اریبی برای هر روش در بودجه‌ی ثابت شات.
- **Definitions:** Pauli-sum observable، shot budget، classical shadows، RMSE، bias.

### 2. Overall Description
- **Product Functions:** تولید/بارگذاری هامیلتونی، اجرای سه تخصیص‌دهنده، محاسبه‌ی خطا نسبت به مقدار دقیق، گزارش مقایسه‌ای.
- **User Characteristics:** پژوهشگر محاسبات کوانتومی با دانش پایه‌ی پایتون و مفاهیم اندازه‌گیری.
- **Constraints:** فقط شبیه‌سازی؛ بدون اتصال به سخت‌افزار کوانتومی واقعی.
- **Assumptions:** هامیلتونی‌ها به‌صورت مجموع جملات پائولی با ضرایب حقیقی نمایش داده می‌شوند.

### 3. Specific Requirements
- **Functional:**
  - FR-01: بارگذاری هامیلتونی از فایل JSON با ساختار `{n_qubits, terms:[{pauli:"XIZ", coeff:0.5}]}`.
  - FR-02: تولید حالت‌های ورودی تصادفی (statevector) با seed قابل تنظیم.
  - FR-03: اجرای تخمین با تخصیص یکنواخت: تقسیم مساوی شات بین جملات.
  - FR-04: اجرای تخمین با تخصیص وزن‌دار: تخصیص متناسب با \(|c_i|\) یا واریانس جمله.
  - FR-05: اجرای تخمین با classical shadows (نمونه‌گیری تصادفی از پایه‌های پائولی).
  - FR-06: محاسبه‌ی RMSE و اریبی برای هر روش و ذخیره در CSV.
  - FR-07: تولید نمودار مقایسه‌ای (اختیاری).
- **External Interfaces:** ورودی JSON؛ خروجی CSV و لاگ متنی.
- **Performance:** اجرای یک آزمایش با ۳ کیوبیت و ۱۲ جمله در کمتر از ۵ ثانیه روی CPU معمولی.
- **Design Constraints:** پایتون ۳.۱۰+، NumPy، SciPy؛ بدون وابستگی به سخت‌افزار کوانتومی.

### 4. قابلیت ثبت اختراع و نوآوری پیشنهادی
- **ارزیابی اولیه:** قابلیت ثبت اختراع **پایین**. تخصیص شات، classical shadows، و مقایسه‌ی آن‌ها در ادبیات فنی به‌طور گسترده پوشش داده شده‌اند. با این حال، یک شکاف بالقوه وجود دارد: **تخصیص تطبیقی بودجه بر اساس عدم قطعیت بیزی پسین و ساختار گروه‌بندی جملات**.
- **نوآوری پیشنهادی برای SRS:** افزودن یک ماژول «تخصیص تطبیقی» که پس از یک دور اندازه‌گیری مقدماتی، واریانس پسین هر جمله را تخمین می‌زند و بودجه‌ی باقیمانده را به‌صورت بهینه (مثلاً با حل یک برنامه‌ریزی محدب کوچک) بین جملات یا گروه‌های جملات تقسیم می‌کند. ادعای روشی بالقوه: «تخصیص پویای شات با استفاده از گراف اشتراک کیوبیت و به‌روزرسانی بیزی، با معیار حداقل‌سازی RMSE مورد انتظار در بودجه‌ی ثابت.»
- **اقدام لازم:** جستجوی دقیق پیشینه در پایگاه‌های پتنت (Google Patents, Espacenet) با کلیدواژه‌های “adaptive shot allocation Pauli grouping Bayesian” و مشورت با وکیل.

---

## پیوست: کد تولید داده مصنوعی

```python
#!/usr/bin/env python3
"""تولید داده مصنوعی برای qshot-lab: هامیلتونی‌های پائولی و حالت‌های تصادفی."""
import json
import numpy as np
import os

def random_pauli_hamiltonian(n_qubits, n_terms, seed=42, max_coeff=1.0):
    rng = np.random.default_rng(seed)
    paulis = ['I', 'X', 'Y', 'Z']
    terms = []
    for _ in range(n_terms):
        s = ''.join(rng.choice(paulis) for _ in range(n_qubits))
        if s == 'I' * n_qubits:
            s = 'Z' + s[1:]  # جلوگیری از جمله همانی خالص
        coeff = float(rng.uniform(-max_coeff, max_coeff))
        if abs(coeff) < 1e-6:
            coeff = 0.1
        terms.append({"pauli": s, "coeff": coeff})
    return {"n_qubits": n_qubits, "terms": terms}

def random_statevector(n_qubits, seed=7):
    rng = np.random.default_rng(seed)
    dim = 2 ** n_qubits
    real = rng.normal(size=dim)
    imag = rng.normal(size=dim)
    vec = real + 1j * imag
    return (vec / np.linalg.norm(vec)).tolist()

def exact_pauli_expectation(state, pauli_string):
    """محاسبه دقیق <ψ|P|ψ> برای رشته پائولی."""
    psi = np.array(state, dtype=complex)
    mats = {
        'I': np.eye(2),
        'X': np.array([[0, 1], [1, 0]]),
        'Y': np.array([[0, -1j], [1j, 0]]),
        'Z': np.array([[1, 0], [0, -1]]),
    }
    op = np.array([[1.0 + 0j]])
    for ch in pauli_string:
        op = np.kron(op, mats[ch])
    return float(np.real(psi.conj() @ op @ psi))

if __name__ == "__main__":
    out_dir = "data/qshot-lab"
    os.makedirs(out_dir, exist_ok=True)
    for n_q, n_t in [(3, 12), (4, 20), (5, 30)]:
        ham = random_pauli_hamiltonian(n_q, n_t, seed=100 + n_q)
        with open(f"{out_dir}/hamiltonian_{n_q}q.json", "w") as f:
            json.dump(ham, f, indent=2)
        st = random_statevector(n_q, seed=200 + n_q)
        with open(f"{out_dir}/state_{n_q}q.json", "w") as f:
            json.dump(st, f)
        rows = []
        for term in ham["terms"]:
            val = exact_pauli_expectation(st, term["pauli"])
            rows.append((term["pauli"], term["coeff"], val, term["coeff"] * val))
        np.savetxt(
            f"{out_dir}/exact_terms_{n_q}q.csv",
            rows,
            delimiter=",",
            header="pauli,coeff,expectation,weighted",
            fmt="%s",
            comments="",
        )
    print(f"داده‌های qshot-lab در {out_dir} تولید شد.")
