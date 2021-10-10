import { MigrationDefinitionWithName } from "umzug/lib/migrationsList";
import { zap } from "../helpers/sql";

const migration: MigrationDefinitionWithName = {
	name: "02-board_sharings",
	async up() {
		await zap.transaction.single/*sql*/`
            create table 
                board_sharings
            ( 
                id uuid primary key,
                "boardId" uuid not null references boards(id),
                "userId" uuid not null references users(id),
                "createdAt" timestamptz not null,
                "updatedAt" timestamptz not null,
				"isArchived" boolean not null
            );
        `;
	},
	async down() {
		await zap.transaction.single/*sql*/`
            drop table board_sharings; 
        `;
	},
};

export default migration;
