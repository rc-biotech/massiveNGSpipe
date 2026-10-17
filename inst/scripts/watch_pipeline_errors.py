#!/usr/bin/env python3
"""Follow massiveNGSpipe error events and print one line per new event (or forward them).

The pipeline appends JSON lines to <project>/error_logs/events.jsonl (R/pipeline_errors.R): event = error,
resolved, crash or stage_gave_up. This watcher tails those files for several projects at once, starting at the
current end (old events are not replayed unless --replay N is given), and prints a compact line per event, so it
can run under a notifier such as Claude Code's Monitor or `watch`. Standard library only.

  python3 watch_pipeline_errors.py                         # all NGS_pipeline* and local_pipeline_* projects
  python3 watch_pipeline_errors.py --projects ~/livemount/Bio_data/NGS_pipeline_RNA --replay 5
  python3 watch_pipeline_errors.py --agent-sync claude-megafst --agent-sync-cli /path/agent_sync.py
  python3 watch_pipeline_errors.py --once                  # print current error flags per project and exit

--agent-sync sends each event as an agent_sync message (kind alert) to that agent as well as printing it.
"""
import argparse, glob, json, os, subprocess, sys, time

DEF_BASE = os.path.expanduser("~/livemount/Bio_data")


def projects_default():
    pats = ["NGS_pipeline", "NGS_pipeline_*", "local_pipeline_*"]
    out = []
    for p in pats:
        out += [d for d in glob.glob(os.path.join(DEF_BASE, p)) if os.path.isdir(d)]
    return sorted(set(out))


def fmt(ev, project_dir):
    proj = ev.get("project") or os.path.basename(project_dir)
    kind = ev.get("event", "?")
    exp, step = ev.get("exp", ""), ev.get("step", "")
    msg = (ev.get("message") or "").replace("\n", " ")[:300]
    hint = ev.get("hint", "")
    tail = f" | hint: {hint}" if hint and hint != "unknown" else ""
    return f"[{ev.get('time', '')}] {proj} {kind.upper()} {exp} ({step}): {msg}{tail}"


def current_flags(project_dir):
    """Error flag files of the newest session (the ones a running pipeline shows)."""
    sessions = sorted(d for d in glob.glob(os.path.join(project_dir, "error_logs", "*")) if os.path.isdir(d))
    if not sessions:
        return None, []
    s = sessions[-1]
    return os.path.basename(s), sorted(os.path.basename(f)[:-4] for f in glob.glob(os.path.join(s, "*.rds")))


def send_agent_sync(cli, agent, line, ev):
    subj = f"pipeline {ev.get('event', '?')}: {ev.get('exp', '')} ({ev.get('step', '')})"[:120]
    try:
        subprocess.run([sys.executable, cli, "send", "--from", "mngsp-error-watch", "--to", agent, "--kind", "alert",
                        "--subject", subj, "--body", line], check=False, timeout=30,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:  # never stop watching because a notification failed
        print(f"(agent_sync send failed: {e})", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--projects", nargs="*", help="project directories (default: all NGS_pipeline*/local_pipeline_*)")
    ap.add_argument("--poll", type=float, default=5.0, help="seconds between checks (default 5)")
    ap.add_argument("--replay", type=int, default=0, help="also print the last N existing events per project")
    ap.add_argument("--once", action="store_true", help="print current error flags per project and exit")
    ap.add_argument("--agent-sync", metavar="AGENT", help="also send each event to this agent_sync agent")
    ap.add_argument("--agent-sync-cli", default=os.path.expanduser(
        "~/shared_workspace/hakon_playground/agent_sync/agent_sync.py"), help="path to agent_sync.py")
    a = ap.parse_args()
    projects = [os.path.expanduser(p) for p in (a.projects or projects_default())]

    if a.once:
        for p in projects:
            sess, flags = current_flags(p)
            if sess:
                print(f"{os.path.basename(p)} session {sess}: {len(flags)} error flag(s)"
                      + (": " + ", ".join(flags) if flags else ""))
        return

    offsets = {}
    for p in projects:
        f = os.path.join(p, "error_logs", "events.jsonl")
        if os.path.exists(f):
            if a.replay > 0:
                with open(f, encoding="utf-8", errors="replace") as fh:
                    lines = fh.readlines()[-a.replay:]
                for l in lines:
                    try:
                        print(fmt(json.loads(l), p), flush=True)
                    except ValueError:
                        pass
            offsets[f] = os.path.getsize(f)
        else:
            offsets[f] = 0
    print(f"watching {len(projects)} project(s) for pipeline error events", flush=True)
    while True:
        for f in list(offsets):
            try:
                size = os.path.getsize(f)
            except OSError:
                continue
            if size < offsets[f]:  # truncated or replaced: start over from the beginning
                offsets[f] = 0
            if size == offsets[f]:
                continue
            with open(f, encoding="utf-8", errors="replace") as fh:
                fh.seek(offsets[f])
                chunk = fh.read()
            # only consume complete lines; a partial last line is read again next round
            end = chunk.rfind("\n")
            if end < 0:
                continue
            offsets[f] += len(chunk[:end + 1].encode("utf-8"))
            for l in chunk[:end + 1].splitlines():
                if not l.strip():
                    continue
                try:
                    ev = json.loads(l)
                except ValueError:
                    continue
                line = fmt(ev, os.path.dirname(os.path.dirname(f)))
                print(line, flush=True)
                if a.agent_sync:
                    send_agent_sync(a.agent_sync_cli, a.agent_sync, line, ev)
        time.sleep(a.poll)


if __name__ == "__main__":
    main()
