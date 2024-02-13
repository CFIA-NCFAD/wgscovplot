import {createEffect, createSignal} from "solid-js";
import {setState} from "../state";

const initializeTheme = () => {
  let theme;
  if (typeof localStorage !== "undefined" && localStorage.getItem("theme")) {
    theme = localStorage.getItem("theme") as "light" | "dark";
  } else if (window.matchMedia("(prefers-color-scheme: dark)").matches) {
    theme = "dark";
  } else {
    theme = "light";
  }
  return theme;
};
export const ThemeButton = () => {
  const [theme, setTheme] = createSignal<string>(initializeTheme());

  createEffect(() => {
    const root = document.documentElement;
    if (theme() === "light") {
      root.classList.remove("dark");
      localStorage.setItem("theme", "light");
    } else {
      root.classList.add("dark");
      localStorage.setItem("theme", "dark");
    }
  });

  return <button class="btn rounded dark:bg-slate-500 bg-slate-700 hover:bg-slate-900 p-1 dark:hover:bg-slate-400"
                 title={theme() === "light" ? "Switch to Dark Mode!" : "Switch to Light Mode!"}
                 type="button" onClick={() => {
    setTheme((t) => (t === "light" ? "dark" : "light"));
    setState("chartOptions", "darkMode", theme() !== "light");
  }}>
    {theme() === "light" ? "🌙" : "☀️"}
  </button>
}