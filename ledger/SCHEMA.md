# Reading the ledger

`campaign.jsonl` is one JSON object per line, appended in the order events happened and never
rewritten. 384 lines, 26 campaigns, 25 September to 4 October 2026.

Fields present on most events:

| field | meaning |
|---|---|
| `type` | the event class, listed below |
| `campaign_id` | the campaign this event belongs to |
| `campaign` | present on the event that opens a campaign; carries `campaign_id` |

Event classes, with the count in this ledger:

| type | n | what it records |
|---|---|---|
| `derived` | 114 | a result computed from a previous one |
| `round` | 96 | one cycle of the inner loop |
| `dispatched` | 79 | a role was given work |
| `status_changed` | 28 | a campaign or theory changed state |
| `created` | 26 | a campaign was opened |
| `believed` | 18 | a belief was declared, carrying theories |
| `objective_falsifier_declared` | 7 | a falsifier was named before the measurement that would answer it |
| `stopping_revised` | 6 | the stopping condition was changed |
| `objective_declared` | 6 | an objective was set |
| `redescribed` | 2 | a campaign was restated |
| `columns_declared` | 1 | an output schema was fixed |
| `judged` | 1 | an adjudication |

Figure 1 counts only `round`, `believed`, `objective_falsifier_declared`, `dispatched` and `judged`,
and only for events naming a campaign that was opened. `verify/verify_ledger.py` applies exactly
that filter and prints the result beside the figure's numbers.

The four excerpts quoted verbatim in Figure 1F are drawn from the agent transcript rather than from
this ledger, except the first, which is a ledger entry.
